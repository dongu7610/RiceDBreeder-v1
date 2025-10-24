#==============
# pedi_module.py — CSV-only minimal pedigree helper (no DB)
import os
import re
import pandas as pd
import networkx as nx
from typing import Tuple, Set, List, Dict, Any

# -------------------------------------------------------
# 1) CSV 경로: base_dir/ordi_arc3/*.csv  (요청 사항 반영)
# -------------------------------------------------------
BASE_DIR   = os.path.dirname(os.path.abspath(__file__))
    
DATA_DIR = os.path.join(BASE_DIR, "data")

CSV_PHENO    = os.path.join(DATA_DIR, "phenotype_data_en.csv")   # rice_characteristics 대체
CSV_RESOURCE = os.path.join(DATA_DIR, "eval_value_list4.csv")    # rice_resource_types 대체
CSV_VCFMAP   = os.path.join(DATA_DIR, "result_nabic2.csv")       # rice_vcf_mapping 대체
CSV_THAL3    = os.path.join(DATA_DIR, "osa_fin1234.csv")         # thal3_data 대체
# 인코딩 기본값(필요시 조정)
ENC_RESOURCE = "cp949"
ENC_VCFMAP   = "cp949"
ENC_THAL3    = "utf-8"


# -------------------------------------------------------
# 2) 유틸
# -------------------------------------------------------
def _read_csv_safe(path: str, **kwargs) -> pd.DataFrame:
    if not os.path.exists(path):
        return pd.DataFrame()
    try:
        return pd.read_csv(path, **kwargs)
    except Exception:
        # 인코딩 문제 대비 등
        return pd.read_csv(path)

def _clean_id(x) -> str:
    s = str(x) if x is not None else ""
    return re.sub(r"\s+", "", s).upper()

def _normalize_no_space(s) -> str:
    if not isinstance(s, str):
        return ""
    return re.sub(r"\s+", "", s.lower())

def _normalize_byeo_drop_no_space(s) -> str:
    if not isinstance(s, str):
        return ""
    s = s.lower().replace("byeo", "")
    return re.sub(r"\s+", "", s)

def _compact_line_name(s) -> str:
    """라인명 숫자+공백 패턴 보정 (예: 'X 123' → 'X123'), 공백제거"""
    if not isinstance(s, str):
        return ""
    s = s.lower()
    s = re.sub(r"(\d+)\s+", r"\1", s)
    return re.sub(r"\s+", "", s)

# -------------------------------------------------------
# 3) CSV 로더: resource / vcf / thal3
# -------------------------------------------------------
def _load_resource() -> pd.DataFrame:
    raw = _read_csv_safe(CSV_RESOURCE, encoding=ENC_RESOURCE, dtype=str)
    if raw.empty:
        return pd.DataFrame(columns=["id", "variety_name_en", "line_name_en"])

    # 컬럼 맵핑 (한국어/영문 혼용 대비)
    colmap = {
        "자원번호": "id",
        "자원명(품종명영문)": "variety_name_en",
        "자원명(계통명영문)": "line_name_en",
        "variety_name_en": "variety_name_en",
        "line_name_en": "line_name_en",
        "id": "id",
    }
    for c in list(colmap.keys()):
        if c in raw.columns:
            raw.rename(columns={c: colmap[c]}, inplace=True)

    for need in ["id", "variety_name_en", "line_name_en"]:
        if need not in raw.columns:
            raw[need] = None

    df = raw[["id", "variety_name_en", "line_name_en"]].copy()
    df["id"] = df["id"].map(_clean_id)
    # 정규화 키 준비 (매칭용)
    df["_var_no_space"]      = df["variety_name_en"].map(_normalize_no_space)
    df["_var_byeo_no_space"] = df["variety_name_en"].map(_normalize_byeo_drop_no_space)
    df["_line_compact"]      = df["line_name_en"].map(_compact_line_name)
    return df.drop_duplicates(subset=["id"])

def _load_vcfmap() -> pd.DataFrame:
    raw = _read_csv_safe(CSV_VCFMAP, encoding=ENC_VCFMAP, dtype=str)
    if raw.empty:
        return pd.DataFrame(columns=["variety_id", "vcf_tag"])

    colmap = {
        "자원번호": "variety_id",
        "Entry_No.": "vcf_tag",
        "variety_id": "variety_id",
        "vcf_tag": "vcf_tag",
    }
    for c in list(colmap.keys()):
        if c in raw.columns:
            raw.rename(columns={c: colmap[c]}, inplace=True)

    for need in ["variety_id", "vcf_tag"]:
        if need not in raw.columns:
            raw[need] = None

    df = raw[["variety_id", "vcf_tag"]].copy()
    df["variety_id"] = df["variety_id"].map(_clean_id)
    df["vcf_tag"] = df["vcf_tag"].astype(str).str.strip()
    df = df[df["variety_id"].ne("")].drop_duplicates("variety_id")
    return df

def _load_thal3_edges() -> pd.DataFrame:
    raw = _read_csv_safe(CSV_THAL3, encoding=ENC_THAL3, dtype=str)
    if raw.empty:
        return pd.DataFrame(columns=["cultivar_name", "parent"])

    # 컬럼 리네임
    colmap = {"Cultivar_Name": "cultivar_name", "P": "parent"}
    for c in list(colmap.keys()):
        if c in raw.columns:
            raw.rename(columns={c: colmap[c]}, inplace=True)
    for need in ["cultivar_name", "parent"]:
        if need not in raw.columns:
            raw[need] = None

    df = raw[["cultivar_name", "parent"]].dropna()
    df = df[(df["cultivar_name"].astype(str) != "") & (df["parent"].astype(str) != "")]
    return df


# -------------------------------------------------------
# 4) Minimal PedigreeApp (CSV only)
#    - 공개 메서드:
#      * get_variety_info_by_pedigree(pedigree_name)
#      * get_connected_nodes(center_node, levels_up, levels_down)
#      * create_cytoscape_elements(nodes_set, edges_set)
# -------------------------------------------------------
class PedigreeApp:
    def __init__(self):
        # CSV 로드
        self.resource_df = _load_resource()                 # id / variety_name_en / line_name_en / 정규화키
        self.vcf_df      = _load_vcfmap()                   # variety_id / vcf_tag
        self.thal3_df    = _load_thal3_edges()              # parent -> cultivar_name

        # VCF 매핑 dict
        self.vcf_map = self.vcf_df.set_index("variety_id")["vcf_tag"].to_dict()

        # NetworkX DiGraph 구성 (parent → child)
        self.graph = nx.DiGraph()
        for _, r in self.thal3_df.iterrows():
            p = str(r["parent"]).strip()
            c = str(r["cultivar_name"]).strip()
            if p and c:
                self.graph.add_edge(p, c)

    # -------------------------
    # 내부 헬퍼: 리소스에서 id/vcf 찾기
    # -------------------------
    def _match_resource_row(self, pedigree_name: str):
        """
        pedigree_name(= variety 이름)으로 resource_df 행을 찾는다.
        매칭 규칙:
          - variety_name_en: no_space / byeo-drop-no-space 비교
          - line_name_en: compact 비교
        """
        if not pedigree_name:
            return None

        key_no_space = _normalize_no_space(pedigree_name)
        key_byeo     = _normalize_byeo_drop_no_space(pedigree_name)
        key_line     = _compact_line_name(pedigree_name)

        df = self.resource_df

        hits = df[
            (df["_var_no_space"] == key_no_space) |
            (df["_var_byeo_no_space"] == key_byeo) |
            (df["_line_compact"] == key_line)
        ]
        if hits.empty:
            return None
        return hits.iloc[0]  # 첫 매치 우선

    # -------------------------
    # 공개 API 1: IT/VCF 조회
    # -------------------------
    def get_variety_info_by_pedigree(self, pedigree_name: str) -> Dict[str, Any]:
        """
        입력된 계보명(pedigree_name)으로부터 IT 번호와 VCF 상태를 CSV에서 조회.
        반환: {"IT Number": <str|None>, "VCF Status": <str|None|'No VCF'>}
        """
        row = self._match_resource_row(pedigree_name)
        if row is None:
            return {"IT Number": None, "VCF Status": None}

        it_number = row["id"]
        vcf = self.vcf_map.get(it_number, None)
        vcf_status = vcf if (vcf is not None and vcf != "") else "No VCF"
        return {"IT Number": it_number, "VCF Status": vcf_status}

    # -------------------------
    # 공개 API 2: 연결 노드 탐색
    # -------------------------
    def get_connected_nodes(self, center_node: str, levels_up: int, levels_down: int) -> Tuple[Set[str], Set[Tuple[str, str]]]:
        """
        center_node 기준으로 위(조상) levels_up 단계, 아래(자손) levels_down 단계 탐색
        반환: (nodes_set, edges_set)
        """
        nodes: Set[str] = set([center_node])
        edges: Set[Tuple[str, str]] = set()

        # 위로(부모) BFS
        from collections import deque
        q = deque([(center_node, 0)])
        while q:
            node, level = q.popleft()
            if level < levels_up:
                for parent in self.graph.predecessors(node):
                    nodes.add(parent)
                    edges.add((parent, node))
                    q.append((parent, level + 1))

        # 아래로(자식) BFS
        q = deque([(center_node, 0)])
        while q:
            node, level = q.popleft()
            if level < levels_down:
                for child in self.graph.successors(node):
                    nodes.add(child)
                    edges.add((node, child))
                    q.append((child, level + 1))

        return nodes, edges

    # (선택) 외부에서 부모/자식 쓰려면 사용
    def _get_node_parents_children(self, node_id: str) -> Tuple[List[str], List[str]]:
        parents = list(self.graph.predecessors(node_id))
        children = list(self.graph.successors(node_id))
        return parents, children

    # -------------------------
    # 공개 API 3: Cytoscape elements 생성
    # -------------------------
    def create_cytoscape_elements(self, nodes_set: Set[str], edges_set: Set[Tuple[str, str]]) -> List[Dict[str, Any]]:
        """
        외부에서 받은 (nodes, edges)로 Cytoscape elements 구성
        """
        elements: List[Dict[str, Any]] = []

        # 노드
        for node in nodes_set:
            varinfo = self.get_variety_info_by_pedigree(str(node))
            it_number = varinfo.get("IT Number")
            vcf_status = varinfo.get("VCF Status")

            has_it = it_number is not None and it_number != "No information"
            has_vcf = (vcf_status is not None) and (vcf_status != "No VCF")
                    # ✅ color_class 결정
            if has_vcf and has_it:
                color_class = "both"
            elif has_vcf:
                color_class = "vcf"
            elif has_it:
                color_class = "it"
            else:
                color_class = "none"

            parents, children = self._get_node_parents_children(node)

            elements.append({
                "data": {
                    "id": str(node),
                    "label": str(node),
                    "it_number": it_number if has_it else "No information",
                    "vcf_status": vcf_status if has_vcf else "No VCF",
                    "has_it": has_it,
                    "has_vcf": has_vcf,
                    "color_class": color_class, 
                    "parents": parents,
                    "children": children,
                },
                "classes": "normal",
            })

        # 에지
        for s, t in edges_set:
            elements.append({"data": {"source": str(s), "target": str(t)}})

        return elements
    

    # pedi_module.py 내 PedigreeApp 클래스 안
    def list_all_cultivars(self):
        """
        CSV(thal3_df)에서 cultivar_name, parent 컬럼을 모두 읽어 전체 목록 반환
        """
        if getattr(self, "thal3_df", None) is None or self.thal3_df.empty:
            return []

        cols = []
        for col in ["cultivar_name", "parent"]:
            if col in self.thal3_df.columns:
                vals = self.thal3_df[col].dropna().astype(str).str.strip()
                cols.extend(vals.tolist())

        return sorted(set([n for n in cols if n]))