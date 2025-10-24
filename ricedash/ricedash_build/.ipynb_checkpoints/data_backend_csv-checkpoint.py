# data_backend_csv.py
# CSV 백엔드: 기존 PostgreSQL 쿼리를 pandas 로 1:1 치환
import os
import re
import pandas as pd
from typing import List, Dict, Optional, Tuple

# ====== 경로 설정 ======
# 필요시 절대경로로 교체하세요
BASE_DIR = os.path.dirname(__file__)


CSV_PHENO    = os.getenv("CSV_PHENO",    os.path.join(BASE_DIR, "data", "phenotype_data_en.csv"))
CSV_RESOURCE = os.getenv("CSV_RESOURCE", os.path.join(BASE_DIR, "data", "eval_value_list4.csv"))
CSV_VCFMAP   = os.getenv("CSV_VCFMAP",   os.path.join(BASE_DIR, "data", "result_nabic2.csv"))
CSV_THAL3    = os.getenv("CSV_THAL3",    os.path.join(BASE_DIR, "data", "osa_fin1234.csv"))


# ====== 로드(전역 캐시) ======
def _read_csv_safe(path, **kwargs):
    if not os.path.exists(path):
        # 존재하지 않으면 빈 DF 반환 (코드가 깨지지 않도록)
        return pd.DataFrame()
    return pd.read_csv(path, **kwargs)

# rice_characteristics 대체: 표현형 테이블
PHENO_DF = _read_csv_safe(CSV_PHENO, dtype=str, encoding="utf-8")

# rice_resource_types 대체
#   - id, resource_type, variety_name_en, line_name_en 로 정규화
_resource_raw = _read_csv_safe(CSV_RESOURCE, encoding="cp949")
def _clean_id(x):
    s = str(x).strip().replace(" ", "").upper()
    return s

if not _resource_raw.empty:
    _resource = pd.DataFrame({
        "id": _resource_raw.get("자원번호", pd.Series(dtype=str)).map(_clean_id),
        "resource_type": _resource_raw.get("자원구분", pd.Series(dtype=str)),
        "variety_name_en": _resource_raw.get("자원명(품종명영문)", pd.Series(dtype=str)).where(lambda s: ~s.isna(), None),
        "line_name_en": _resource_raw.get("자원명(계통명영문)", pd.Series(dtype=str)).where(lambda s: ~s.isna(), None)
    })
else:
    _resource = pd.DataFrame(columns=["id","resource_type","variety_name_en","line_name_en"])

# rice_vcf_mapping 대체: variety_id <-> vcf_tag(Entry_No)
_vcfmap_raw = _read_csv_safe(CSV_VCFMAP, encoding="cp949") 
if not _vcfmap_raw.empty:
    _vcfmap = pd.DataFrame({
        "variety_id": _vcfmap_raw.get("자원번호", pd.Series(dtype=str)).map(_clean_id),
        "vcf_tag": _vcfmap_raw.get("Entry_No.", pd.Series(dtype=str)).astype(str).str.strip()
    }).dropna(subset=["variety_id"]).drop_duplicates("variety_id")
else:
    _vcfmap = pd.DataFrame(columns=["variety_id","vcf_tag"])

# thal3_data 대체: pedigree 이름 매칭용
_thal3_raw = _read_csv_safe(CSV_THAL3 ,encoding="utf-8")
if not _thal3_raw.empty:
    # 컬럼 리네임
    _thal3 = _thal3_raw.rename(columns={
        "Cultivar_Name": "cultivar_name",
        "P": "parent"
    })
    for c in ("cultivar_name","parent"):
        if c not in _thal3.columns:
            _thal3[c] = None
    _thal3["cultivar_name_l"] = _thal3["cultivar_name"].astype(str).str.lower().str.replace(r"\s+","",regex=True)
    _thal3["parent_l"]        = _thal3["parent"].astype(str).str.lower().str.replace(r"\s+","",regex=True)
else:
    _thal3 = pd.DataFrame(columns=["cultivar_name","parent","cultivar_name_l","parent_l"])

# ====== Category/Type 매핑 (DB 없이 하드코딩) ======
def get_all_category_mappings() -> Dict[str, Dict]:
    return {
        'yield_quality': {
            'name': 'yield_quality',
            'traits': [
                'endosperm_type', 'hull_color', 'thousand_grain_weight', 'amylose',
                'protein', 'alkali_digestion', 'grain_length', 'grain_width',
                'grain_shape', 'white_core', 'white_belly', 'grains_per_panicle',
                'amylose_nir', 'protein_nir', 'total_phenol_ug_per_g_gae', 'total_anthocyanin',
                'abt_scavenging_ug_per_ml', 'dpph_scavenging'
            ]
        },
        'cultivation': {
            'name': 'cultivation',
            'traits': [
                'planting_date', 'seeding_date', 'heading_date', 
                'variety_type', 'panicles_per_plant'
            ]
        },
        'disaster_resistance': {
            'name': 'disaster_resistance',
            'traits': [
                'pi39_resistance', 'cold_resistance', 'bacterial_leaf_blight_k1',
                'stripe_disease', 'brown_planthopper', 'pit_resistance',
                'bacterial_leaf_blight_k3a', 'pita_pita2_resistance', 'pid_t2_resistance',
                'low_temp_germination', 'low_temp_germination_grade', 'viviparous',
                'pita_resistance', 'pii_resistance', 'bakanae_disease', 'pib_resistance',
                'bacterial_leaf_blight_k3', 'pi5_resistance', 'pik_resistance',
                'pikp_resistance', 'pizt_resistance', 'bacterial_leaf_blight',
                'bacterial_leaf_blight_k2', 'leaf_blast', 'pikm_resistance',
                'shattering_resistance'
            ]
        },
        'morphological': {
            'name': 'morphological',
            'traits': [
                'leaf_color', 'tillering_habit', 'internode_color', 'flag_leaf_angle',
                'panicle_exertion', 'culm_length', 'panicle_length', 'awn_length',
                'awn_color', 'leaf_sheath_color', 'leaf_erectness', 'apiculus_color',
                'apiculus_color_mature'
            ]
        }
    }

def determine_trait_type_fallback(trait_name: str) -> str:
    date_traits = ['planting_date', 'seeding_date', 'heading_date']
    if trait_name in date_traits:
        return 'date'
    numeric_traits = ['low_temp_germination','amylose_nir','abt_scavenging_ug_per_ml',
                      'total_phenol_ug_per_g_gae','dpph_scavenging','amylose',
                      'total_anthocyanin','thousand_grain_weight','protein','protein_nir',
                      'culm_length','panicle_length']
    if trait_name in numeric_traits:
        return 'numeric'
    return 'categorical'

# ====== DB 대체 함수들 ======

def get_resource_types() -> List[str]:
    if _resource.empty:
        return []
    # 번역 테이블(미상→Unknown 등)이 필요하면 여기서 처리
    mapping = {
        'nan': 'Unknown', '미상': 'Unknown', '유전재료': 'Germplasm',
        '육성계통': 'Breeding Line', '육성품종': 'Cultivar',
        '잡초형': 'Weedy Type', '재래종': 'Landrace'
    }
    rt = _resource['resource_type'].astype(str).map(lambda v: mapping.get(v, v))
    return sorted(rt.dropna().unique().tolist())

def get_db_data(resource_types: Optional[List[str]] = None) -> pd.DataFrame:
    """
    기존:
      SELECT rc.* FROM rice_characteristics rc
      WHERE rc.id IN (SELECT id FROM rice_resource_types WHERE resource_type IN ...)
    대체: PHENO_DF + _resource 조인 후 resource_type 필터
    """
    if PHENO_DF.empty:
        return pd.DataFrame()
    ph = PHENO_DF.copy()
    if 'id' not in ph.columns:
        # id 컬럼 없으면 실패 방지
        return pd.DataFrame()

    # 정규화: id 키
    ph['id'] = ph['id'].map(_clean_id)

    # 리소스 타입 정규화
    res = _resource.copy()
    mapping = {
        'nan': 'Unknown', '미상': 'Unknown', '유전재료': 'Germplasm',
        '육성계통': 'Breeding Line', '육성품종': 'Cultivar',
        '잡초형': 'Weedy Type', '재래종': 'Landrace'
    }
    res['resource_type'] = res['resource_type'].astype(str).map(lambda v: mapping.get(v, v))

    if resource_types:
        res = res[res['resource_type'].isin(resource_types)]

    # JOIN
    df = ph.merge(res[['id','resource_type','variety_name_en','line_name_en']],
                  on='id', how='inner')

    # variety_type 컬럼 필요 시 채움
    if 'variety_type' not in df.columns:
        df['variety_type'] = df['resource_type']
    return df

def get_available_columns(df: pd.DataFrame) -> List[Dict[str,str]]:
    """
    기존 SQL(category_mappings + categories JOIN) 제거.
    - 카테고리: get_all_category_mappings()에서 찾고
    - 없으면 'Other'
    """
    excluded = {'id', 'variety_type'}
    catmap = get_all_category_mappings()
    trait2cat = {}
    for cat, blob in catmap.items():
        for t in blob['traits']:
            trait2cat[t] = cat

    cols = []
    for col in df.columns:
        if col in excluded: 
            continue
        label = f"{trait2cat.get(col, 'Other')}: {col}"
        cols.append({"label": label, "value": col})
    return cols

def get_column_type_from_db(column_name: str) -> str:
    """
    DB 없이 타입 판정. (fallback 규칙)
    """
    return determine_trait_type_fallback(column_name)

# ====== 결과 테이블(복잡한 SQL) 대체 ======

def _normalize_name(s: Optional[str]) -> Optional[str]:
    if s is None or pd.isna(s): 
        return None
    s = str(s).strip()
    if not s:
        return None
    return s

def _match_thal3(variety_name_en: Optional[str], line_name_en: Optional[str]) -> Optional[str]:
    """
    기존 SQL 서브쿼리들을 판다스로 재현:
    - cultivar_name / parent 컬럼에 대해 lower/공백제거/특정 치환(byeo 등)하여 매칭
    """
    if _thal3.empty:
        return None

    def _norm(v: Optional[str], drop_byeo=False) -> str:
        if v is None: 
            return ""
        s = str(v).lower()
        if drop_byeo:
            s = s.replace("byeo", "")
        s = re.sub(r"\s+", "", s)
        return s

    # 1) variety_name_en 우선
    if variety_name_en:
        key = _norm(variety_name_en, drop_byeo=False)
        key_nospace = key
        key_dropbyeo = _norm(variety_name_en, drop_byeo=True)

        # cultivar_name == key
        hit = _thal3.loc[_thal3["cultivar_name_l"] == key]
        if not hit.empty:
            return hit["cultivar_name"].iloc[0]

        # cultivar_name == nospace
        hit = _thal3.loc[_thal3["cultivar_name_l"] == key_nospace]
        if not hit.empty:
            return hit["cultivar_name"].iloc[0]

        # cultivar_name == dropbyeo
        hit = _thal3.loc[_thal3["cultivar_name_l"] == key_dropbyeo]
        if not hit.empty:
            return hit["cultivar_name"].iloc[0]

        # parent == key / nospace / dropbyeo
        for col in ["parent_l"]:
            for k in (key, key_nospace, key_dropbyeo):
                hit = _thal3.loc[_thal3[col] == k]
                if not hit.empty:
                    # parent 일치 시 parent 반환
                    return _thal3["parent"].iloc[hit.index[0]]

    # 2) line_name_en 보정: 숫자+공백 정리
    if line_name_en:
        line_norm = re.sub(r"(\d+)\s+", r"\1", str(line_name_en).lower()).strip()
        line_norm = re.sub(r"\s+","", line_norm)

        hit = _thal3.loc[_thal3["cultivar_name_l"] == line_norm]
        if not hit.empty:
            return hit["cultivar_name"].iloc[0]
        hit = _thal3.loc[_thal3["parent_l"] == line_norm]
        if not hit.empty:
            return _thal3["parent"].iloc[0]

    return None

def build_results_table_data(filtered_ids: List[str]) -> pd.DataFrame:
    """
    기존 거대 SQL을 판다스로 재현.
    출력 컬럼:
      IT Number, Rice Variety Name, Rice Line Name, Variety_type, VCF Status,
      processed_name, detection_status
    """
    if not filtered_ids:
        return pd.DataFrame(columns=[
            "IT Number","Rice Variety Name","Rice Line Name",
            "Variety_type","VCF Status","processed_name","detection_status"
        ])

    # JOIN: resource + pheno + vcf
    df = get_db_data(None)  # 전체
    df = df[df["id"].isin(filtered_ids)].copy()

    # Suwon -> Suweon 교정
    df["line_name_en"] = df["line_name_en"].apply(
        lambda x: re.sub(r"^Suwon", "Suweon", x) if isinstance(x, str) else x
    )

    # VCF 상태 조인
    vmap = _vcfmap.set_index("variety_id")["vcf_tag"]
    df["VCF Status"] = df["id"].map(vmap).fillna("No VCF")

    # processed_name 계산 (thal3 매칭)
    df["processed_name"] = [
        _match_thal3(_normalize_name(vn), _normalize_name(ln))
        for vn, ln in zip(df.get("variety_name_en"), df.get("line_name_en"))
    ]
    df["processed_name"] = df["processed_name"].fillna("nopedi")

    # detection_status
    def _det(s, vn, ln):
        if s and s != "nopedi":
            if vn or ln:
                return "Pedigree"
        return "Not found"
    df["detection_status"] = [
        _det(pn, vn, ln) for pn, vn, ln in zip(df["processed_name"], df.get("variety_name_en"), df.get("line_name_en"))
    ]

    out = pd.DataFrame({
        "IT Number": df["id"],
        "Rice Variety Name": df.get("variety_name_en"),
        "Rice Line Name": df.get("line_name_en"),
        "Variety_type": df.get("variety_type"),
        "VCF Status": df["VCF Status"].replace("No VCF", None),
        "processed_name": df["processed_name"],
        "detection_status": df["detection_status"]
    }).sort_values("IT Number")

    return out
