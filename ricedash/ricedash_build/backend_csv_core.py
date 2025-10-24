# -*- coding: utf-8 -*-
"""
CSV 백엔드 코어: DB 없이 CSV에서 로드하여 판다스로 처리.
- load_resource_types():  id, variety_name_en, line_name_en
- load_vcf_mapping():     variety_id, vcf_tag
- load_thal3():           cultivar_name, parent, feature (+ 정규화 컬럼)
- (선택) load_phenotype_info(): phenotype_data_en.csv 그대로 로드
"""

import os
import re
import pandas as pd
from typing import Optional

# ===== 경로 설정 =====
BASE_DIR     = os.path.dirname(os.path.abspath(__file__))
CSV_PHENO    = os.getenv("CSV_PHENO",    os.path.join(BASE_DIR, "data", "phenotype_data_en.csv"))
CSV_RESOURCE = os.getenv("CSV_RESOURCE", os.path.join(BASE_DIR, "data", "eval_value_list4.csv"))
CSV_VCFMAP   = os.getenv("CSV_VCFMAP",   os.path.join(BASE_DIR, "data", "result_nabic2.csv"))
CSV_THAL3    = os.getenv("CSV_THAL3",    os.path.join(BASE_DIR, "data", "osa_fin1234.csv"))

# ===== 인코딩 기본 =====
ENC_PHENO    = os.getenv("ENC_PHENO",    "utf-8")
ENC_RESOURCE = os.getenv("ENC_RESOURCE", "cp949")
ENC_VCFMAP   = os.getenv("ENC_VCFMAP",   "cp949")
ENC_THAL3    = os.getenv("ENC_THAL3",    "utf-8")

# ===== 공용 유틸 =====
def _read_csv_safe(path: str, encoding: Optional[str] = None, **kwargs) -> pd.DataFrame:
    """
    CSV 읽기: 없으면 빈 DF, 인코딩 실패 시 기본(dtypes=str)로 재시도
    """
    if not os.path.exists(path):
        return pd.DataFrame()
    try:
        return pd.read_csv(path, encoding=encoding, dtype=str, **kwargs)
    except Exception:
        # encoding 지정 실패 시 마지막 방어
        return pd.read_csv(path, dtype=str, **kwargs)

def _clean_id(x) -> str:
    """
    ID 정규화: 문자열화 + 공백 제거 + 대문자
    """
    if x is None:
        return ""
    try:
        if pd.isna(x):
            return ""
    except Exception:
        pass
    s = str(x).strip()
    s = re.sub(r"\s+", "", s)
    return s.upper()

def _safe_str(x) -> str:
    """
    어떤 값이 와도 문자열로 안전 변환 (NaN/None → "")
    """
    if x is None:
        return ""
    try:
        if pd.isna(x):
            return ""
    except Exception:
        pass
    return str(x)

# ===== 로더들 =====
def load_phenotype_info() -> pd.DataFrame:
    """
    rice_characteristics 대체: phenotype_data_en.csv 그대로 반환 (dtype=str)
    """
    df = _read_csv_safe(CSV_PHENO, encoding=ENC_PHENO)
    return df

def load_resource_types() -> pd.DataFrame:
    """
    rice_resource_types 대체 CSV → [id, variety_name_en, line_name_en]
    - 원본 컬럼:
       자원번호 → id
       자원명(품종명영문) → variety_name_en
       자원명(계통명영문) → line_name_en
    """
    df = _read_csv_safe(CSV_RESOURCE, encoding=ENC_RESOURCE)
    if df.empty:
        return pd.DataFrame(columns=["id", "variety_name_en", "line_name_en"])

    colmap = {
        "자원번호": "id",
        "자원명(품종명영문)": "variety_name_en",
        "자원명(계통명영문)": "line_name_en",
        # 이미 영문 컬럼명으로 되어있는 경우를 대비
        "id": "id",
        "variety_name_en": "variety_name_en",
        "line_name_en": "line_name_en",
    }
    for src, dst in colmap.items():
        if src in df.columns:
            df.rename(columns={src: dst}, inplace=True)

    # 필수 컬럼 보장
    for k in ["id", "variety_name_en", "line_name_en"]:
        if k not in df.columns:
            df[k] = None

    # 문자열/NaN 안전화
    df["id"] = df["id"].map(_clean_id)
    df["variety_name_en"] = df["variety_name_en"].map(_safe_str)
    df["line_name_en"]    = df["line_name_en"].map(_safe_str)

    return df[["id", "variety_name_en", "line_name_en"]].drop_duplicates()

def load_vcf_mapping() -> pd.DataFrame:
    """
    rice_vcf_mapping 대체 CSV → [variety_id, vcf_tag]
    - 원본 컬럼:
      자원번호 → variety_id
      Entry_No. → vcf_tag
    """
    df = _read_csv_safe(CSV_VCFMAP, encoding=ENC_VCFMAP)
    if df.empty:
        return pd.DataFrame(columns=["variety_id", "vcf_tag"])

    colmap = {
        "자원번호": "variety_id",
        "Entry_No.": "vcf_tag",
        # 이미 영문 컬럼명인 경우
        "variety_id": "variety_id",
        "vcf_tag": "vcf_tag",
    }
    for src, dst in colmap.items():
        if src in df.columns:
            df.rename(columns={src: dst}, inplace=True)

    # 필수 컬럼 보장 & 정규화
    for k in ["variety_id", "vcf_tag"]:
        if k not in df.columns:
            df[k] = None

    df["variety_id"] = df["variety_id"].map(_clean_id)
    df["vcf_tag"]    = df["vcf_tag"].map(_safe_str).str.strip()

    # 유효한 variety_id만
    df = df[(df["variety_id"] != "") & (df["variety_id"] != "-")]

    # 한 id당 하나의 vcf_tag만 남김
    return df[["variety_id", "vcf_tag"]].drop_duplicates("variety_id")

def load_thal3() -> pd.DataFrame:
    """
    thal3_data 대체 CSV → [cultivar_name, parent, feature, cultivar_name_l, parent_l]
    - 정규화 컬럼(cultivar_name_l, parent_l)은 소문자+공백제거
    """
    df = _read_csv_safe(CSV_THAL3, encoding=ENC_THAL3)
    if df.empty:
        return pd.DataFrame(columns=["cultivar_name", "parent", "feature",
                                     "cultivar_name_l", "parent_l"])

    # 컬럼 리네임 (원본 헤더 다양성 대비)
    colmap = {
        "Cultivar_Name": "cultivar_name",
        "cultivar_name": "cultivar_name",
        "P": "parent",
        "parent": "parent",
        "feature": "feature",
    }
    for src, dst in colmap.items():
        if src in df.columns:
            df.rename(columns={src: dst}, inplace=True)

    for k in ["cultivar_name", "parent", "feature"]:
        if k not in df.columns:
            df[k] = None

    df["cultivar_name"] = df["cultivar_name"].map(_safe_str)
    df["parent"]        = df["parent"].map(_safe_str)
    df["feature"]       = df["feature"].map(_safe_str)

    # 정규화 열 추가 (소문자 + 공백제거)
    def _norm_no_space(s: str) -> str:
        s = _safe_str(s).lower()
        return re.sub(r"\s+", "", s)

    df["cultivar_name_l"] = df["cultivar_name"].map(_norm_no_space)
    df["parent_l"]        = df["parent"].map(_norm_no_space)

    return df[["cultivar_name", "parent", "feature", "cultivar_name_l", "parent_l"]].drop_duplicates()
