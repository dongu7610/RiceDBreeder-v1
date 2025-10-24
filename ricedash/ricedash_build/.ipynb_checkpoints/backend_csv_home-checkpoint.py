# -*- coding: utf-8 -*-
"""
홈(app_home) 전용 헬퍼:
- build_variety_df(): 홈 드롭다운/카드용 DF 생성
  (id, variety_name_en, line_name_en, display_name, vcf_status, has_pedigree, pedigree_name)
- 내부적으로 load_resource_types / load_vcf_mapping / load_thal3 사용
- 문자열/NaN 안전화 처리 포함 (replace 호출 에러 방지)
"""

import re
import pandas as pd
from typing import Optional

from backend_csv_core import (
    load_resource_types,
    load_vcf_mapping,
    load_thal3,
    _safe_str,
)

# ---------- 정규화/매칭 유틸 ----------

def _normalize_text_no_space(s: str) -> str:
    """소문자 + 공백 제거. s가 문자열이 아니어도 안전."""
    s = _safe_str(s).lower()
    return re.sub(r"\s+", "", s)

def _compact_line_name(s: str) -> str:
    """라인명 숫자+공백 패턴 보정 (예: 'X 123' → 'X123'), 소문자 + 전체 공백 제거"""
    s = _safe_str(s).lower()
    s = re.sub(r"(\d+)\s+", r"\1", s)
    return re.sub(r"\s+", "", s)

def _find_pedigree_name(variety_name_en: Optional[str],
                        line_name_en: Optional[str],
                        thal_df: pd.DataFrame) -> Optional[str]:
    """
    DB의 서브쿼리 매칭을 pandas 로 치환:
    - cultivar_name / parent 에 대해:
        1) variety_name_en 기반: 공백삭제키, byeo 제거 후 공백삭제키
        2) line_name_en 기반: 숫자-공백 보정키
    - 순서대로 히트하면 해당 원본명(cultivar_name 또는 parent) 반환
    """
    if thal_df is None or thal_df.empty:
        return None

    # 정규화 컬럼이 없다면 생성(안전)
    t = thal_df.copy()
    if "cultivar_name_l" not in t.columns:
        t["cultivar_name_l"] = t["cultivar_name"].map(_normalize_text_no_space)
    if "parent_l" not in t.columns:
        t["parent_l"] = t["parent"].map(_normalize_text_no_space)

    # variety 기반 키
    v = _safe_str(variety_name_en)
    v_key0 = _normalize_text_no_space(v)                       # 공백삭제
    v_key1 = _normalize_text_no_space(v.replace("byeo", ""))   # 'byeo' 제거 + 공백삭제

    # line 기반 키
    l = _safe_str(line_name_en)
    l_key  = _compact_line_name(l)

    # 탐색 헬퍼
    def _try_hit(keys):
        for k in keys:
            if not k:
                continue
            # cultivar_name 우선
            hit = t.loc[t["cultivar_name_l"] == k]
            if not hit.empty:
                return hit["cultivar_name"].iloc[0]
            # parent 대체
            hit = t.loc[t["parent_l"] == k]
            if not hit.empty:
                return hit["parent"].iloc[0]
        return None

    # 1) variety_name 기반
    cand = _try_hit([v_key0, v_key1])
    if cand:
        return cand

    # 2) line_name 기반
    cand = _try_hit([l_key])
    if cand:
        return cand

    return None

# ---------- 홈용 DF 빌더 ----------

def build_variety_df() -> pd.DataFrame:
    """
    홈(검색 카드)에서 쓰는 데이터프레임 구성:
      id, variety_name_en, line_name_en, display_name, vcf_status, has_pedigree, pedigree_name
    """
    rt = load_resource_types()   # id / variety_name_en / line_name_en
    vm = load_vcf_mapping()      # variety_id / vcf_tag
    th = load_thal3()            # cultivar_name / parent / feature (+ normalized)

    if rt.empty:
        return pd.DataFrame(columns=[
            "id", "variety_name_en", "line_name_en",
            "vcf_status", "has_pedigree", "pedigree_name", "display_name"
        ])

    # join VCF
    merged = rt.merge(vm, how="left", left_on="id", right_on="variety_id")
    merged["vcf_status"] = merged["vcf_tag"].where(merged["vcf_tag"].notna(), "No VCF")

    # 계보 이름 탐색 (NaN 안전)
    merged["pedigree_name"] = merged.apply(
        lambda r: _find_pedigree_name(r.get("variety_name_en"), r.get("line_name_en"), th),
        axis=1
    )
    merged["has_pedigree"] = merged["pedigree_name"].notna()

    # display_name
    def _disp(r):
        idv = r.get("id")
        vn  = _safe_str(r.get("variety_name_en"))
        ln  = _safe_str(r.get("line_name_en"))
        if vn and ln:
            return f"{idv} - {vn} / {ln}"
        if vn:
            return f"{idv} - {vn}"
        if ln:
            return f"{idv} - {ln}"
        return _safe_str(idv)

    merged["display_name"] = merged.apply(_disp, axis=1)

    out = merged[[
        "id", "variety_name_en", "line_name_en",
        "vcf_status", "has_pedigree", "pedigree_name", "display_name"
    ]].drop_duplicates()

    return out.sort_values("id").reset_index(drop=True)
