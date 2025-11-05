import os
import sys
import json
import time
from urllib.parse import unquote
from typing import Union, Iterable
import copy
import networkx as nx
import dash
import dash_bootstrap_components as dbc
from dash import html, dcc, Input, Output, State, callback_context, no_update, callback, ALL, dash_table,ctx,MATCH,clientside_callback
from dash.exceptions import PreventUpdate
import dash_cytoscape as cyto
import pandas as pd
import psycopg2
from psycopg2 import pool
import plotly.graph_objects as go
import numpy as np
import plotly.express as px
import re
from datetime import datetime
from typing import Dict, List, Optional, Tuple
from functools import lru_cache
from backend_csv_core import load_resource_types, load_vcf_mapping

from multi_filter import (
    FilterManager, get_initial_filter_values,
     apply_stored_filter
)
try:
    from backend_csv_core import load_phenotype_info  # 존재하지 않으면 except로 넘어감
    _HAS_PHENO = True
except Exception:
    _HAS_PHENO = False

_BASE_DIR = os.path.dirname(os.path.abspath(__file__))
from data_backend_csv import (
    get_db_data, 
    get_available_columns,
    get_resource_types,   
    build_results_table_data,    
    get_column_type_from_db,
)
from backend_csv_home import build_variety_df
from backend_csv_core import (
    load_resource_types,
    load_vcf_mapping,
    load_thal3,
    _safe_str,
)


def valid_vcf_value(val):
    if not val:
        return None
    val = str(val).strip()
    if val.lower() in {"none", "no vcf", "n/a", "null", "no information"}:
        return None
    return val


def _legend_item(color, shape, label, border="#333"):
    """범례용 아이콘 (색상/모양 선택 가능)"""
    return html.Div(
        [
            html.Div(
                style={
                    "width": "14px",
                    "height": "14px",
                    "backgroundColor": color,
                    "border": f"1.5px solid {border}",
                    "borderRadius": "50%" if shape == "circle" else "3px",
                    "display": "inline-block",
                    "marginRight": "3px",
                }
            ),
            html.Span(label, style={"fontSize": "10px", "color": "#f8f9fa"}),
        ],
        style={"display": "flex", "alignItems": "center", "marginRight": "10px"},
    )

def _legend_ring_item(color, shape, label):
    """중앙이 비어 있고 외곽선만 있는 ring형 아이콘"""
    return html.Div(
        [
            html.Div(
                style={
                    "width": "14px",
                    "height": "14px",
                    "backgroundColor": "transparent",
                    "border": f"2.2px solid {color}",
                    "borderRadius": "50%" if shape == "circle" else "3px",
                    "display": "inline-block",
                    "marginRight": "3px",
                    "boxShadow": f"0 0 3px {color}",  # 약간의 glow 효과
                }
            ),
            html.Span(label, style={"fontSize": "10px", "color": "#f8f9fa"}),
        ],
        style={"display": "flex", "alignItems": "center", "marginRight": "10px"},
    )



def get_color_class(has_vcf: bool, has_it: bool) -> str:
    """
    노드의 has_vcf / has_it 상태를 기반으로 color_class를 반환한다.
    반환값: "both" | "vcf" | "it" | "none"
    """
    if has_vcf and has_it:
        return "both"
    elif has_vcf:
        return "vcf"
    elif has_it:
        return "it"
    else:
        return "none"

pill_css = {
    "display":"flex","gap":"6px","flexWrap":"wrap","alignItems":"center",
    "border":"1px dashed #dfe6ec","padding":"8px","borderRadius":"8px","background":"#fbfdff"
}

base_table_props = dict(
    row_selectable='multi',     # ✅ 좌측 체크박스
    selected_rows=[],           # 초기 선택 없음
    editable=False,             # 데이터 편집 불필요
    cell_selectable=False,
    style_cell={'textAlign':'left','padding':'8px','fontSize':'12px'},
    style_header={'backgroundColor':'rgb(235,235,235)','fontWeight':'bold'},
    page_size=5,
    css=[                       # ✅ 좌측 선택 체크박스 열 스타일(폭/정렬)
        {"selector": ".dash-select-cell", "rule": "width: 42px; text-align: center;"},
        {"selector": ".dash-spreadsheet-container .dash-spreadsheet-inner table", "rule":"width:100%;"},
    ]
)

chip_style = {
    "display":"inline-flex","alignItems":"center","gap":"6px",
    "padding":"4px 8px","borderRadius":"999px",
    "background":"#f1f3f5","border":"1px solid #dee2e6",
    "fontSize":"12px","lineHeight":"18px"
}
chip_x_style = {
    "cursor":"pointer","border":"none","background":"transparent",
    "fontWeight":"bold","fontSize":"14px","color":"#c92a2a"
}

def is_nopedi_case(fixed_vcf, fixed_vcf_nopedi):
    """현재 상태가 'nopedi'인지 판별"""
    if fixed_vcf_nopedi and not fixed_vcf:
        return True
    return False



def make_alert_callbacks(store_id, span_id, timer_id, label_text, color):
    """
    공통 Alert 콜백 생성기
    - store_id: 감시할 store
    - span_id: 메시지 표시할 span ID
    - timer_id: interval ID
    - label_text: 메시지 내용
    - color: 강조 색상
    """

    # --- 1️⃣ Alert 표시 ---
    @callback(
        Output(span_id, 'children'),
        Output(span_id, 'style'),
        Output(timer_id, 'disabled'),
        Input(store_id, 'data'),
        prevent_initial_call=True
    )
    def show_alert(data):
        if not data:
            return "", {'opacity': '0'}, True

        # snp-occurrence-store 전용 — enabled=False이면 알람 X
        if store_id == 'snp-occurrence-store':
            enabled = data.get("enabled", False)
            threshold = data.get("threshold", None)
            if not enabled or threshold is None:
                # 초기 렌더나 비활성화 상태에서는 알람 표시 안 함
                raise dash.exceptions.PreventUpdate

        return (
            f"{label_text}",
            {
                'opacity': '1',
                'color': color,
                'fontWeight': 'bold',
                'fontSize': '12px',
                'transition': 'opacity 0.5s ease'
            },
            False  # 타이머 활성화
        )

    # --- 2️⃣ Alert 숨김 ---
    @callback(
        Output(span_id, 'style', allow_duplicate=True),
        Output(timer_id, 'disabled', allow_duplicate=True),
        Input(timer_id, 'n_intervals'),
        State(timer_id, 'disabled'),
        prevent_initial_call='initial_duplicate'
    )
    def hide_alert(n_intervals, timer_disabled):
        if timer_disabled:
            raise dash.exceptions.PreventUpdate
        return (
            {'opacity': '0', 'transition': 'opacity 0.5s ease'},
            True
        )

def get_default_stylesheet2(basenode: str = None):
    """
    기본 Cytoscape 스타일시트.
    has_vcf / has_it 조합에 따라 색상을 다르게 적용하고,
    basenode(기준 노드)가 있으면 추가 강조 스타일을 부여한다.
    """
    base_styles = [
        # -------------------------------------------------
        # 🔹 모든 노드의 기본 스타일
        # -------------------------------------------------
        {
            'selector': 'node',
            'style': {
                'content': 'data(label)',
                'text-valign': 'center',
                'text-halign': 'center',
                #'background-color': '#6FB1FC',  # 기본값 (fallback)
                'width': 70,
                'height': 70,
                'font-size': 13,
                'fontWeight': 'bold',
                'font-family': 'Arial, sans-serif',
                'color': '#FFFFFF',
                'text-outline-color': '#000000',
                'text-outline-width': 1.5,
                'border-width': 3,
                'border-color': '#2C3E50',
                'text-wrap': 'wrap',
                'text-max-width': 68,
                'transition-property': 'background-color, border-color',
                'transition-duration': '0.5s',
            },
        },
        # -------------------------------------------------
        # 🔸 에지 기본 스타일
        # -------------------------------------------------
        {
            'selector': 'edge',
            'style': {
                'curve-style': 'bezier',
                'target-arrow-shape': 'triangle',
                'arrow-scale': 1.4,
                'line-color': '#ccc',
                'target-arrow-color': '#ccc',
            },
        },

        # -------------------------------------------------
        # 🎨 color_class별 색상 정의
        # -------------------------------------------------
        # 🔵 VCF only
        {
            'selector': 'node[color_class = "vcf"]',
            'style': {
                'background-color': '#1C6BA0',  # 파랑
                'border-color': '#2c3e50',
            }
        },
        # 🟢 IT only
        {
            'selector': 'node[color_class = "it"]',
            'style': {
                'background-color': '#27ae60',  # 초록
                'border-color': '#1e8449',
            }
        },
        # 🟠 VCF + IT both
        {
            'selector': 'node[color_class = "both"]',
            'style': {
                'background-color': '#e67e22',  # 주황
                'border-color': '#a35400',
            }
        },
        # ⚪ none (데이터 없음)
        {
            'selector': 'node[color_class = "none"]',
            'style': {
                'background-color': '#bdc3c7',  # 회색
                'border-color': '#7f8c8d',
            }
        },
    ]

    # -------------------------------------------------
    # ⭐ 기준 노드 강조 (Shadow + Underlay)
    # -------------------------------------------------
    if basenode:
        base_styles.append({
            'selector': f'node[id = "{basenode}"]',
            'style': {
                "border-width": 5,
                "border-color": "#000000",
                "shadow-color": "#7f8c8d",
                "shadow-blur": 12,
                "shadow-opacity": 0.6,
                "shadow-offset-x": 3,
                "shadow-offset-y": 3,
                "underlay-color": "#ecf0f1",
                "underlay-opacity": 0.9,
                "underlay-padding": 4,
            }
        })

    return base_styles


def _legend_circle(fill_color, border_color, label_text):
    """작은 원 + 텍스트 조합 생성기"""
    return html.Div([
        html.Div(style={
            'width': '12px',
            'height': '12px',
            'borderRadius': '50%',
            'backgroundColor': fill_color,
            'display': 'inline-block',
            'marginRight': '6px',
            'border': f'1px solid {border_color}'
        }),
        html.Span(label_text, style={'fontSize': '11px'})
    ], style={'display': 'flex', 'alignItems': 'center', 'marginBottom': '4px'})

def build_stylesheet_with_available2(
    *,
    base=None,
    available_ids=None,
    expanded_ids=None,
    expanded_child_ids=None,
    base_name=None,
    add_group_labels=None
):
    """
    계층형 스타일 통합 (v3):
    - available: 빨강 내부
    - base/add: 짙은 회색 외곽, (suffix 포함)
    - base/add + available: 내부 빨강 + 외곽 짙은 회색 + suffix 라벨
    - expand: 초록 테두리
    - expand-child: 주황 내부
    - expand+child: 초록 외곽 + 주황 내부
    """
    base = base or []
    ss = get_default_stylesheet2()

    available_ids = set(available_ids or [])
    expanded_ids = set(expanded_ids or [])
    expanded_child_ids = set(expanded_child_ids or [])

    dark_gray = "#424949"
    darker_border = "#212121"

    # ✅ helper: content 구성
    def _label_content(name: str, suffix: str | None = None):
        return f"{name}\n({suffix})" if suffix else name

    # --- 1️⃣ 기본 available (빨강 내부)
    for nid in available_ids:
        ss.append({
            "selector": f'node[id = "{_esc_attr_value(nid)}"]',
            "style": {
                "background-color": "#e74c3c",
                "border-color": "#c0392b",
                "border-width": 3,
                "color": "#fff",
            },
        })

    # --- 2️⃣ expand / expand-child (주황 + 초록)
    for nid in expanded_child_ids | expanded_ids:
        if nid in available_ids:
            continue
        has_expand = nid in expanded_ids
        has_child = nid in expanded_child_ids
        style = {}
        if has_expand and has_child:
            style.update({
                "background-color": "#f39c12",
                "border-color": "#27ae60",
                "border-width": 4,
            })
        elif has_expand:
            style.update({
                "border-color": "#27ae60",
                "border-width": 4,
            })
        elif has_child:
            style.update({
                "background-color": "#f39c12",
                "border-color": "#e67e22",
                "border-width": 3,
            })
        if style:
            ss.append({"selector": f'node[id = "{_esc_attr_value(nid)}"]', "style": style})

    # --- 3️⃣ Base Node 처리
    if base_name:
        style = {
            "color": "#fff",
            "font-weight": "600",
            "font-size": 12,
            "text-wrap": "wrap",
            "text-max-width": 68,
            "content": _label_content(base_name, "default"),
        }
        if base_name in available_ids:
            style.update({
                "background-color": "#e74c3c",  # 내부 빨강
                "border-color": darker_border,   # 외곽 짙은 회색
                "border-width": 4,
            })
        else:
            style.update({
                "background-color": dark_gray,
                "border-color": darker_border,
                "border-width": 3,
            })
        ss.append({"selector": f'node[id = "{_esc_attr_value(base_name)}"]', "style": style})

    # --- 4️⃣ Add Group 처리
    if add_group_labels:
        for name, suffix in add_group_labels.items():
            if name == base_name:
                continue
            style = {
                "color": "#fff",
                "font-weight": "500",
                "font-size": 12,
                "text-wrap": "wrap",
                "text-max-width": 68,
                "content": _label_content(name, suffix),
            }
            if name in available_ids:
                style.update({
                    "background-color": "#e74c3c",   # 내부 빨강
                    "border-color": darker_border,   # 외곽 짙은 회색
                    "border-width": 4,
                })
            else:
                style.update({
                    "background-color": dark_gray,   # 내부 회색
                    "border-color": darker_border,
                    "border-width": 3,
                })
            ss.append({"selector": f'node[id = "{_esc_attr_value(name)}"]', "style": style})

    return ss + base






MAX_SEL = 200
# ---- create_filter_component_fix.py (또는 당신의 페이지 파일 상단) ----

def _pedigree_to_bool(x):
    if pd.isna(x): return False
    s = str(x).strip().lower()
    # build_results_table_data()의 'Pedigree' / 'Not found' 기준 대응
    return s in ("pedigree","true","1","y","yes","found","detected")

def _is_vcf_ok(val: str) -> bool:
    if val is None:
        return False
    s = str(val).strip()
    if s == "":
        return False
    # "no vcf", "nan", "none", "null" → False
    if re.fullmatch(r"(?i)\s*(no\s*vcf|nan|none|null)\s*", s):
        return False
    return True

def create_filter_component_fix(
    filter_id_num: int,
    column: str,
    column_type: str,
    df: pd.DataFrame,
    id_suffix: str = "_fix",
):
    """
    기존 컴포넌트 구조는 유지하되, 패턴 매칭 ID의 'type' 필드에만 suffix를 부여하여
    동일 앱 내 다른 페이지/섹션과 ID 충돌을 방지합니다.
    사용 예) new_filter = create_filter_component_fix(n, col, col_type, df, id_suffix="_fix")
    """
    preview_data = df[column].head(10).tolist() if column in df.columns else []  # 필요 시 사용
    months = [{"label": f"{i}", "value": i} for i in range(1, 13)]
    days = [{"label": f"{i}", "value": i} for i in range(1, 32)]
    init = get_initial_filter_values(df, column_type, column)

    # numeric 파라미터
    if column_type == "numeric":
        rng = init.get("range", [0, 100])
        rmin, rmax = float(rng[0]), float(rng[1])
        rstep = (rmax - rmin) / 100 if rmin != rmax else 1
    else:
        rmin = 0
        rmax = 100
        rstep = 1

    t = lambda base: f"{base}{id_suffix}"  # 'type' 값에 suffix 부여

    label = f"Filter : {column}"
    return html.Div(
        [
            html.Div(
                [
                    html.H4(
                        label,
                        style={
                            "display": "inline-block",
                            "fontSize": "16px",
                            "color": "#2c3e50",
                            "fontWeight": "bold",
                            "margin": "0",
                        },
                    ),
                    html.Button(
                        "×",
                        id={"type": t("remove-filter"), "index": filter_id_num},
                        style={
                            "float": "right",
                            "border": "none",
                            "background": "none",
                            "fontSize": "20px",
                            "cursor": "pointer",
                            "color": "#dc3545",
                            "fontWeight": "bold",
                        },
                    ),
                ],
                style={"marginBottom": "15px"},
            ),

            # DATE
            html.Div(
                [
                    html.Div(
                        [
                            html.Label(
                                "Start date:",
                                style={
                                    "fontWeight": "bold",
                                    "marginRight": "20px",
                                    "width": "80px",
                                    "display": "inline-block",
                                    "verticalAlign": "top",
                                    "paddingTop": "5px",
                                    "fontSize": "14px",
                                    "color": "#495057",
                                },
                            ),
                            html.Div(
                                [
                                    html.Div(
                                        [
                                            html.Label("Month:"),
                                            dcc.Dropdown(
                                                id={
                                                    "type": t("filter-startdate-month"),
                                                    "index": filter_id_num,
                                                },
                                                options=months,
                                                value=init.get("start_month"),
                                                placeholder="month",
                                                style={"width": "100px"},
                                            ),
                                        ],
                                        style={
                                            "display": "inline-block",
                                            "marginRight": "10px",
                                        },
                                    ),
                                    html.Div(
                                        [
                                            html.Label("Day:"),
                                            dcc.Dropdown(
                                                id={
                                                    "type": t("filter-startdate-day"),
                                                    "index": filter_id_num,
                                                },
                                                options=days,
                                                value=init.get("start_day"),
                                                placeholder="day",
                                                style={"width": "100px"},
                                            ),
                                        ],
                                        style={"display": "inline-block"},
                                    ),
                                ]
                            ),
                        ],
                        style={"marginBottom": "10px"},
                    ),
                    html.Div(
                        [
                            html.Label(
                                "End date:",
                                style={
                                    "fontWeight": "bold",
                                    "marginRight": "20px",
                                    "width": "80px",
                                    "display": "inline-block",
                                    "verticalAlign": "top",
                                    "paddingTop": "5px",
                                    "fontSize": "14px",
                                    "color": "#495057",
                                },
                            ),
                            html.Div(
                                [
                                    html.Div(
                                        [
                                            html.Label("Month:"),
                                            dcc.Dropdown(
                                                id={
                                                    "type": t("filter-enddate-month"),
                                                    "index": filter_id_num,
                                                },
                                                options=months,
                                                value=init.get("end_month"),
                                                placeholder="month",
                                                style={"width": "100px"},
                                            ),
                                        ],
                                        style={
                                            "display": "inline-block",
                                            "marginRight": "10px",
                                        },
                                    ),
                                    html.Div(
                                        [
                                            html.Label("Day:"),
                                            dcc.Dropdown(
                                                id={
                                                    "type": t("filter-enddate-day"),
                                                    "index": filter_id_num,
                                                },
                                                options=days,
                                                value=init.get("end_day"),
                                                placeholder="day",
                                                style={"width": "100px"},
                                            ),
                                        ],
                                        style={"display": "inline-block"},
                                    ),
                                ]
                            ),
                        ]
                    ),
                ],
                id={"type": t("date-container"), "index": filter_id_num},
                style={"display": "block" if column_type == "date" else "none"},
            ),

            # NUMERIC
            html.Div(
                [
                    dcc.RangeSlider(
                        id={"type": t("numeric-range"), "index": filter_id_num},
                        min=rmin,
                        max=rmax,
                        step=rstep,
                        value=[rmin, rmax],
                        tooltip={"placement": "bottom", "always_visible": True},
                        marks={rmin: f"{rmin:.2f}", rmax: f"{rmax:.2f}"},
                    )
                ],
                id={"type": t("numeric-container"), "index": filter_id_num},
                style={"display": "block" if column_type == "numeric" else "none"},
            ),

            # CATEGORICAL
            html.Div(
                [
                    dcc.Dropdown(
                        id={"type": t("categorical-filter"), "index": filter_id_num},
                        options=[
                            {"label": str(v), "value": str(v)}
                            for v in sorted(df[column].dropna().unique())
                        ]
                        if column in df.columns
                        else [],
                        multi=True,
                        placeholder="Select...",
                        value=(
                            [str(v) for v in sorted(df[column].dropna().unique())]
                            if column in df.columns
                            else []
                        ),
                    )
                ],
                id={"type": t("categorical-container"), "index": filter_id_num},
                style={"display": "block" if column_type == "categorical" else "none"},
            ),

            html.Div(id={"type": t("filter-result"), "index": filter_id_num}),
            html.Hr(),
        ],
        className="filter-block",
    )

def _looks_like_it(s: str) -> bool:
    import re
    return bool(re.match(r"^IT\d+$", s.strip(), re.I))

def _looks_like_entry(s: str) -> bool:
    import re
    # RWG-001, IRIS-xxxx 등 보편 패턴 (필요 시 확장)
    return bool(re.match(r"^[A-Za-z]{2,5}-?\d+$", s.strip()))


def _as_list(v):
    if v is None: return []
    if isinstance(v, (list, tuple, set)): return list(v)
    return [v]

def _uniq(seq):
    seen=set(); out=[]
    for x in seq:
        if x not in seen:
            seen.add(x); out.append(x)
    return out

def load_vcf_info():
    """VCFinfo.csv 파일을 로드하여 DataFrame으로 반환"""
    try:
        BASE_DIR = os.path.dirname(os.path.abspath(__file__))
        vcf_info_path = os.path.join(BASE_DIR, "data", "VCFinfo.csv")
        df = pd.read_csv(vcf_info_path)
        return df
    except Exception as e:
        print(f"Error loading VCFinfo.csv: {e}")
        return pd.DataFrame()

# VCF 정보 전역 로드
VCF_INFO_DF = load_vcf_info()

VCF_MAP_DF  = load_vcf_mapping() 


if not VCF_MAP_DF.empty and "vcf_tag" in VCF_MAP_DF.columns:
    VCF_MAP_DF = VCF_MAP_DF.rename(columns={"vcf_tag": "Entry_No."})

# 하나로 합치기 (Entry_No. 기준)
VCF_ALL_DF = (pd.merge(VCF_INFO_DF, VCF_MAP_DF, on="Entry_No.", how="left")
              if not VCF_MAP_DF.empty else VCF_INFO_DF.copy())

VARIETY_BY_EN_NAME = None
VARIETY_DF: pd.DataFrame = build_variety_df()
def _ensure_variety_by_en_name():
    """영문 품종명 → VARIETY_DF 행 매핑 dict 준비"""
    global VARIETY_BY_EN_NAME
    if VARIETY_BY_EN_NAME is not None:
        return VARIETY_BY_EN_NAME
    df = VARIETY_DF.copy()
    keycol = None
    # display_name에 영문명이 들어있다면 그대로 활용, 없으면 별도 영문 컬럼이 있는지 확인
    for c in ["display_name", "variety_name_en", "Variety Name in English"]:
        if c in df.columns:
            keycol = c; break
    if keycol is None:
        VARIETY_BY_EN_NAME = {}
        return VARIETY_BY_EN_NAME
    m = {}
    for _, r in df.iterrows():
        name = str(r.get(keycol, "")).strip()
        if name:
            m[name.lower()] = r.to_dict()
    VARIETY_BY_EN_NAME = m
    return VARIETY_BY_EN_NAME

def map_vcf_to_it(entry_no: str, en_name: str):
    """
    VCF_INFO_DF 한 행 → (has_pedigree, it_id or None, vcf_status or None)
    - 영문 품종명으로 VARIETY_DF를 역매칭
    """
    en = (en_name or "").strip().lower()
    if not en:
        return False, None, None
    m = _ensure_variety_by_en_name()
    row = m.get(en)
    if not row:
        return False, None, None
    it_id = row.get("id")
    vcf_status = row.get("vcf_status")
    has_pedi = bool(row.get("has_pedigree"))
    return has_pedi, it_id, vcf_status

def attach_vcf_status_to_ids(it_ids):
    """
    IT id 리스트에 vcf_status 붙여서 반환용 dict로 변환
    """
    if not it_ids: return []
    df = VARIETY_DF[VARIETY_DF["id"].isin(it_ids)]
    out=[]
    for _, r in df.iterrows():
        out.append({
            "id": r["id"],
            "display_name": r.get("display_name"),
            "vcf_status": r.get("vcf_status"),
            "has_pedigree": bool(r.get("has_pedigree")),
            "pedigree_name": r.get("pedigree_name") if pd.notna(r.get("pedigree_name")) else None,
        })
    return out

base_table_props2 = dict(
page_size=8,
style_table={"overflowX": "auto"},
style_cell={"fontSize": 12, "whiteSpace": "nowrap", "textOverflow": "ellipsis"},
style_header={"fontWeight": "bold"},
)

_GWAS_STORE = {
    "per_sample": {},   # {sample_id: {status, variant_only(DF), variant_pvalue(DF), error}}
    "combos": {}        # {"S1|S3": {status, requested_samples, matched_samples,
                        #             variant_only(DF), variant_pvalue(DF), error}}
}

LIMIT_MAX = 200  #
opt2_filter_manager = FilterManager()

CATEGORY_MAPPINGS = {
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

selection_stores = html.Div(
[
dcc.Store(id="addv-input-multi", data=[]),
dcc.Store(id="addv-additional-phenotype", data=[]),
dcc.Store(id="addv-additional-vcf", data=[]),
 dcc.Store(id="mannequin-active-filter", data=None),
 dcc.Store(id="opt3-active-filter", data=None),
  dcc.Store(id='vcf-info-global', data=VCF_INFO_DF.to_dict('records')),
   dcc.Store(id="opt3-selected-vcf", data=[]),  
  dcc.Store(id='filter-chain-store_fix', data={"filters": [], "results": []}),
dcc.Store(id="opt2-selected-ids", data=[]),      # persistent selected ITs (IDs)
dcc.Store(id="opt2-last-table-data", data=[]),   # last rendered table rows for mapping id -> cols


dcc.Store(id="addv-input-opt1", data=[]),     # opt1 선택 결과 (IT id list)
dcc.Store(id="addv-input-opt2", data=[]),     # opt2 선택 결과 (IT id list) ← opt2의 addv-input-multi 미러
dcc.Store(id="addv-input-opt3", data=[]),     # opt3 선택 결과 (IT id list)
dcc.Store(id="addv-union", data=[]),          # opt1/2/3 전체 합집합 (UI 요약용)
dcc.Store(id="addv-meta-cache", data=[]),     # union에 대해 조립된 메타 (name/pedigree_label/vcf_status)
# (Apply 결과 저장용 – 선택사항이지만 아래 콜백에 맞춰 추가 권장)
dcc.Store(id="addv-apply-groups-pedi", data={}),
dcc.Store(id="addv-apply-groups-nopedi", data={}),
dcc.Store(id="addv-additional-pedigree", data=[]),
dcc.Store(id="opt3-selected-preview-rows", data=[]),


]
)
def _section_title(text):
    return html.Div([
    html.Div(text, className="text-sm font-semibold", style={"color": "#2c3e50"}),
    html.Hr(style={"margin": "3px 0 12px", "opacity": 0.25}),
    ])
# -----------------------------------------------------------------------------
# Component: Right-side summary panel (sticky inside modal)
# - Input-multiple current values 미리보기
# - Additional(phenotype/vcf) 카운트와 프리뷰
# -----------------------------------------------------------------------------
summary_panel = html.Div(
    [
        _section_title("Current Selection"),
        html.Div(
            [
                html.Div(
                    id="addv-input-multi-preview",
                    style={
                        "maxHeight": "160px",
                        "overflowY": "auto",
                        "border": "1px solid #e1e5ea",
                        "padding": "8px",
                        "borderRadius": "8px",
                        "background": "#fafcfe"
                    }
                ),
                html.Div("Additional lists", className="text-xs",
                         style={"color": "#7f8c8d", "marginTop": "10px"}),

                html.Details([
                    html.Summary("Preview additional lists", className="text-xs",
                                 style={"cursor": "pointer"}),
                    html.Div([
                        html.Div("Phenotype", className="text-xs",
                                 style={"color": "#2c3e50", "marginTop": "6px"}),
                        html.Div(id="addv-additional-phenotype-preview",
                                 style={
                                     "maxHeight": "120px",
                                     "overflowY": "auto",
                                     "border": "1px dashed #e1e5ea",
                                     "padding": "6px",
                                     "borderRadius": "8px"
                                 }),
                        html.Div("VCF", className="text-xs",
                                 style={"color": "#2c3e50", "marginTop": "10px"}),
                        html.Div(id="addv-additional-vcf-preview",
                                 style={
                                     "maxHeight": "120px",
                                     "overflowY": "auto",
                                     "border": "1px dashed #e1e5ea",
                                     "padding": "6px",
                                     "borderRadius": "8px"
                                 }),
                    ])
                ]),

                html.Div(
                    [
                        dbc.Badge("Phenotype: 0", id="addv-additional-phenotype-count",
                                  color="secondary", pill=True, className="me-1"),
                        dbc.Badge("VCF: 0", id="addv-additional-vcf-count",
                                  color="secondary", pill=True),
                        dbc.Badge("Pedigree: 0", id="addv-additional-pedigree-count",
                                  color="secondary", pill=True),
                    ],
                    style={
                        "display": "flex",
                        "gap": "6px",
                        "flexWrap": "wrap",
                        "marginTop": "6px"
                    },
                )
            ]
        ),
    ],
    style={
        "position": "sticky",
        "top": 0,
        "border": "1px solid #e1e5ea",
        "borderLeft": "3px solid #cfd6dc",  # ✅ 구분용 세로선
        "marginLeft": "3px",                # ✅ 좌우 여백
        "marginRight": "3px",
        "borderRadius": "12px",
        "padding": "14px",
        "background": "#ffffff",
        "boxShadow": "0 1px 2px rgba(16,24,40,0.04)"
    }
)
COMMON_OPT_STYLE = {
    "marginTop": "8px",
    "minHeight": "180px",
    "paddingBottom": "10px"
}

# --- Option 1: 버튼 제거 ---
opt1_layout = html.Div(
    [
        _section_title("Search varieties (Pedigree-enabled)"),

        html.Div([
            html.Label("Search varieties", className="text-xs", style={"color": "#2c3e50"}),

            # ✅ opt2의 Quick Add 느낌을 그대로 적용
            html.Div([
                dcc.Dropdown(
                    id="opt1-pedigree-variety-search",
                    options=[],
                    value=None,
                    multi=True,
                    placeholder="Type to search varieties…",
                    maxHeight=280,
                    style={"fontSize": "13px", "minHeight": "40px"},
                ),
            ], style={
                "border": "1px solid #e9ecef",
                "borderRadius": "8px",
                "padding": "6px 8px",
                "background": "#fff",
                "marginTop": "6px",
                "marginBottom": "4px"
            }),

            html.Small(
                "Search or select multiple varieties; results will appear immediately below.",
                style={"color": "#7f8c8d"}
            ),
            html.Div(id="opt1-summary", className="text-xs",
                     style={"color": "#34495e", "marginTop": "8px"})
        ], style=COMMON_OPT_STYLE)
    ]
)


opt2_layout = html.Div(
    [
        _section_title("Phenotype-based checklist (filters)"),

        # Quick input (above resource type)
        html.Div([
                html.Div("Quick Add to Current Selection", className="filter-title"),
                html.Div(id="opt2-quick-panel", style={
                    "display": "flex", "flexWrap": "wrap", "gap": "8px",
                    "minHeight": "32px", "alignItems": "center",
                    "padding": "6px 8px", "border": "1px solid #e9ecef", "borderRadius": "8px",
                    "background": "#fff"
                }),
                html.Small(
                    "Checked items appear here immediately. Up to 200 items are allowed; beyond that you’ll see a warning.",
                    style={"color":"#7f8c8d", "display":"block", "marginTop":"4px"},
                ),
                html.Div(id="opt2-notice", style={"marginTop":"6px"}),  # 경고/안내 표시
            ], className="filter-block", style={"marginBottom":"10px"}),

        html.Div([
            html.Div(
                id="opt2-checklist-count",
                className="text-xs",
                style={"color":"#2c3e50","marginTop":"12px","marginBottom":"6px"}
            ),
            dash_table.DataTable(
                id="opt2-pedi-checklist_fix",
                columns=[
                    {"id":"id","name":"IT Number"},
                    {"id":"name","name":"Variety Name"},
                    {"id":"vcf_status","name":"VCF Status"},
                    {"id":"pedigree_label","name":"Pedigree"},   # ✅ 표시는 이 한 컬럼
                ],
                # hidden_columns 제거 ✅
                data=[],
                row_selectable="multi",
                selected_row_ids=[],
                # row_id 제거 ✅
                **base_table_props2,
            ),
            html.Small(
                "When checked, items are added immediately to ‘Current Selection’.",
                style={"color":"#7f8c8d"}
            ),
        ], style=COMMON_OPT_STYLE),

        # Resource Type Selection (slightly smaller)
        html.Div([
            html.Div("Resource Type Selection", className="filter-title", style={"fontSize":"12px"}),
            dcc.Dropdown(
                id='filter-x_fix',
                options=[{'label': rt, 'value': rt} for rt in get_resource_types()],
                placeholder="Select resource types…",
                multi=True
            ),
            html.Div(id='filter-x-result_fix', className="filter-content"),
        ], className="filter-block", style={"marginBottom":"8px"}),

        # Add filter(+) (shown only when resource selected)
        html.Div([
            dbc.Button(
                "Add filter(+)",
                id="add-filter-btn_fix",
                color="primary",
                size="sm",
                className="add-filter-btn",
                style={'display': 'none'}
            ),
        ], className="filter-btn-container"),

        # Dynamic filters grid
        html.Div(
            id="opt2-phenotype-filter-builder_fix",
            style={
                "border":"1px dashed #dfe6ec",
                "padding":"10px",
                "borderRadius":"8px",
                "background":"#fbfdff",
                "display": "grid",
                "gridTemplateColumns": "repeat(2, minmax(340px, 1fr))",
                "gap": "12px",
            }
        ),

        # Checklist block with count
        

        # Filter picker modal
        dbc.Modal([
            dbc.ModalHeader("Choose Filter"),
            dbc.ModalBody([
                dcc.Dropdown(id='filter-type-selector_fix', placeholder="Select a filter to add")
            ]),
            dbc.ModalFooter([
                dbc.Button("Add", id="confirm-filter-btn_fix", className="ms-auto", color="primary"),
                dbc.Button("Cancel", id="close-modal-btn_fix", className="ms-2")
            ]),
        ], id="filter-modal_fix"),

        # Dynamic filter container (kept for compatibility if you append children here)
        html.Div(id='dynamic-filters_fix'),
    ]
)
# -----------------------------------------------------------------------------
# OPTION 3: VCF group ingestion + checklist handoff
# - VCF 테이블(혹은 업로드/그룹) → 체크리스트
# - 체크된 것은 input-multiple 로, 나머지는 additional-vcf 로
# -----------------------------------------------------------------------------
opt3_layout = html.Div(
    [
        _section_title("Genotype group (pills → checklist)"),

        # (1) VCF 샘플 선택 드롭다운
        html.Div([
                    html.Label("Select VCF samples", className="text-xs"),

                    # ✅ Dropdown 감싸는 border block
                    html.Div([
                        dcc.Dropdown(
                            id="opt3-vcf-sample-dropdown",
                            options=[],
                            placeholder="Choose VCF Entry_No.…",
                            multi=True,
                            clearable=True,
                            style={"fontSize": "13px", "minHeight": "40px"},
                        ),
                    ], style={
                        "border": "1px solid #e9ecef",
                        "borderRadius": "8px",
                        "padding": "6px 8px",
                        "background": "#fff",
                        "marginTop": "6px",
                        "marginBottom": "4px"
                    }),

                    html.Small(
                        "Select one or more VCF samples; these will appear below and sync with the checklist.",
                        style={"color": "#7f8c8d"}
                    ),
                ], style=COMMON_OPT_STYLE),

        # (2) Pill 행 (드롭다운 아래로 이동)
        html.Div(
            [
                html.Div([
                    html.Div("Ecotype", className="text-xs", style={"color":"#2c3e50"}),
                    html.Div(id="ecotype-filter-container2",
                             style={"display":"flex","gap":"6px","flexWrap":"wrap",
                                    "border":"1px dashed #dfe6ec","padding":"8px",
                                    "borderRadius":"8px","background":"#fbfdff"}),
                ], style={"flex":"1 1 280px","minWidth":"260px"}),

                html.Div([
                    html.Div("Variety Group", className="text-xs", style={"color":"#2c3e50"}),
                    html.Div(id="variety-group-filter-container2",
                             style={"display":"flex","gap":"6px","flexWrap":"wrap",
                                    "border":"1px dashed #dfe6ec","padding":"8px",
                                    "borderRadius":"8px","background":"#fbfdff"}),
                ], style={"flex":"2 1 380px","minWidth":"320px"}),
            ],
            style={"display":"flex","gap":"10px","flexWrap":"wrap","marginBottom":"8px"}
        ),

        # (3) 체크리스트 (Pill에 의해 생성 — 체크하면 Dropdown value에 합류)
        html.Div([
            html.Div("Checklist (by Ecotype / Variety Group)", className="text-xs",
                     style={"color":"#2c3e50","marginBottom":"6px"}),
            dash_table.DataTable(
                id="opt3-vcf-checklist",
                columns=[
                    {"id":"Entry_No.","name":"Entry_No."},
                    {"id":"Variety Name in English","name":"Variety (EN)"},
                    {"id":"Ecotype","name":"Ecotype"},
                    {"id":"Variety Group","name":"Variety Group"},
                    {"id":"variety_id","name":"Variety ID"},  # 있으면 표시
                ],
                data=[],
                row_selectable="multi",
                selected_rows=[],
                **base_table_props2,
            ),
            html.Small("When checked, the corresponding Entry_No. will be added to the ‘Select VCF samples’ dropdown.", style={"color":"#7f8c8d"}),
        ], style={"marginBottom":"10px"}),

        # (4) 선택된 VCF 프리뷰 (드롭다운 value 그대로 표시)
        html.Div([
            html.Div("Preview of selected VCF samples", className="text-xs",
                     style={"color":"#2c3e50","marginTop":"6px"}),
            dash_table.DataTable(
                id="opt3-preview-table",
                columns=[
                    {"id":"Entry_No.","name":"Entry_No."},
                    {"id":"Variety Name in English","name":"Variety (EN)"},
                    {"id":"Variety Group","name":"Variety Group"},
                    {"id":"Ecotype","name":"Ecotype"},
                    {"id":"Origin","name":"Origin"},
                    {"id":"NCBI Biosample ID","name":"NCBI Biosample ID"},
                ],
                data=[],
                **base_table_props2,
            ),
        ])
    ]
)



# -----------------------------------------------------------------------------
# MAIN MODAL (single modal with Option1~3)
# - 한 모달 안에서 탭으로 옵션 전환
# -----------------------------------------------------------------------------
additional_varieties_modal = html.Div(
                            [
                            selection_stores,


                            # 트리거 버튼
                            dbc.Button(
                            [html.I(className="fas fa-plus me-1"), "Add Additional Varieties"],
                            id="open-additional-modal-btn",
                            color="primary",
                            outline=False,
                            size="sm",
                            style={"marginBottom": "10px"},
                            ),


                            dbc.Modal(
                            [
                            dbc.ModalHeader(
                            dbc.ModalTitle("Select Additional Varieties for Analysis")
                            ),
                            dbc.ModalBody(
                            html.Div(
                            [
                            # 좌측: 탭, 우측: 요약(sticky)
                            html.Div(
                            [
                            dbc.Tabs(
                            [
                            dbc.Tab(opt1_layout, label="Pedigree Search", tab_id="opt1"),
                            dbc.Tab(opt2_layout, label="Phenotype Filters", tab_id="opt2"),
                            dbc.Tab(opt3_layout, label="Genotype Group", tab_id="opt3"),
                            ],
                            id="addv-tabs",
                            ),
                            ],
                            style={"flex": "3 1 600px", "minWidth": "320px"},
                            ),
                            html.Div(summary_panel, style={"flex": "1 1 280px", "minWidth": "260px"}),
                            ],
                            style={"display": "flex", "gap": "14px", "alignItems": "flex-start", "flexWrap": "wrap",
                                            "maxHeight": "80vh",   # 화면 높이의 80%
                                            "overflowY": "auto",   # 내용이 넘치면 스크롤
                                            "padding": "1rem",},
                            )
                            ),
                            dbc.ModalFooter(
                            [
                            html.Div(id="addv-footer-hint", className="text-xs", style={"color": "#7f8c8d", "marginRight": "auto"}),
                            dbc.Button("Apply", id="addv-apply-btn", color="primary"),
                            dbc.Button("Close", id="addv-close-btn", color="secondary", outline=True),
                            ]
                            ),
                            ],
                            id="additional-varieties-modal",
                            size="xl",
                            is_open=False,
                            scrollable=True,
                            backdrop="static",
                            keyboard=False,
                            centered=True,  
                            #style={"maxWidth": "1100px"},
                            ),
                            ]
)



def _norm_id(x) -> str:
    return str(x).strip().replace(" ", "").upper() if x is not None else ""

def _dash_if_empty(s) -> str:
    s = "" if s is None else str(s).strip()
    return s if s else "-"

def _extract_node_it_map(elements):
    """
    elements에서 '노드'만 고르고 (id, it_number) 맵을 만든다.
    - 노드: data에 id가 있고, source/target이 없음
    - it_number 키는 it_number / IT Number / itNumber 모두 대응
    - 유효한 IT만 수집
    """
    node_ids, it_numbers, node_to_itn = [], [], {}
    for el in elements or []:
        data = el.get("data", {})
        if not isinstance(data, dict):
            continue

        # 노드 판별: id는 있고, edge 키(source/target)는 없어야 함
        is_node = ("id" in data) and ("source" not in data and "target" not in data)
        if not is_node:
            continue

        nid = str(data.get("id")) if data.get("id") is not None else None
        itn = data.get("it_number") or data.get("IT Number") or data.get("itNumber")
        
        itn = _norm_id(itn)

        # 필요시 has_it 플래그까지 체크하려면 주석 해제
        # if data.get("has_it") is not True:
        #     itn = None

        if nid and itn:
            node_ids.append(nid)
            it_numbers.append(itn)
            node_to_itn[nid] = itn

    return node_ids, it_numbers, node_to_itn



def format_group_label(group_type, group_num):
    if group_type == "trait":
        offset = 0
        prefix = "Trait Group"
    elif group_type == "sample":
        offset = 5
        prefix = "Sample Group"
    elif group_type == "variant":
        offset = 10
        prefix = "Variant Group"
    else:
        offset = 0
        prefix = "Group"
    display_num = group_num - offset
    return f"{prefix} {display_num}"

def determine_trait_type_fallback(trait_name):
    """DB 조회 실패시 사용하는 fallback 분류 로직"""
    # 날짜 관련 trait들
    date_traits = ['planting_date', 'seeding_date', 'heading_date']
    if trait_name in date_traits:
        return 'date'
    
    # 정량형으로 알려진 trait들
    numeric_traits = ['low_temp_germination','amylose_nir','abt_scavenging_ug_per_ml','total_phenol_ug_per_g_gae','dpph_scavenging','amylose',
'total_anthocyanin','thousand_grain_weight','protein','protein_nir','culm_length','panicle_length']
    if trait_name in numeric_traits:
        return 'numeric'
    
    # 나머지는 categorical
    return 'categorical'

def get_all_category_mappings():
    """기존 코드와의 호환을 위한 getter"""
    return CATEGORY_MAPPINGS

def get_category_traits(category_name: str):
    """카테고리별 trait 리스트 편의 함수"""
    cm = CATEGORY_MAPPINGS.get(category_name, {})
    return cm.get('traits', [])

def _safe_unique_str(seq):
    seen = set()
    out = []
    for x in seq or []:
        sx = str(x)
        if sx not in seen:
            seen.add(sx)
            out.append(sx)
    return out

def _collect_selected_varieties_data(analysis_nodes):
    """
    analysis_nodes: IT Number 리스트(문자열)
    반환: { node_id(=IT): row(dict) }
    """
    result = {}
    if not analysis_nodes:
        return result
    # 한 번에 슬라이싱
    df = pheno_by_ids(analysis_nodes)
    if df.empty:
        return result
    df = df.astype(str)
    by_id = {str(r['id']): r.to_dict() for _, r in df.iterrows() if 'id' in r}
    for it in analysis_nodes:
        row = by_id.get(str(it))
        if row:
            result[str(it)] = row
    return result


def _pick_rows_by_ids(pheno_df: pd.DataFrame, ids: list[str]) -> dict[str, pd.Series]:
    """
    ids는 IT number(=PHENO_DF['id'])라고 가정.
    존재하지 않는 id는 건너뜁니다.
    """
    out = {}
    if pheno_df is None or pheno_df.empty or not ids:
        return out

    # 문자열 기반 비교 안전화
    df = pheno_df.copy()
    
    if 'id' not in df.columns:
        return out
    df['id'] = df['id'].astype(str) #변경 전: df['id'] == str(node_id)

    for node_id in ids:
        s = df.loc[df['id'] == str(node_id)] #변경 전: df['id'] == str(node_id)
        if not s.empty:
            out[str(node_id)] = s.iloc[0]
    
    return out

def _available_traits(pheno_df: pd.DataFrame, traits: list[str], selected_rows: dict[str, pd.Series]) -> list[str]:
    """
    선택된 품종 중 하나라도 유효한 값이 있는 trait만 남깁니다.
    """
    if pheno_df is None or pheno_df.empty or not traits or not selected_rows:
        return []

    df_cols = set(pheno_df.columns)
    kept = []
    for t in traits:
        if t not in df_cols:
            continue
        has_val = False
        for _, row in selected_rows.items():
            if t in row.index:
                val = row[t]
                if pd.notna(val) and str(val).strip() not in ('', 'nan', 'None', '-'):
                    has_val = True
                    break
        if has_val:
            kept.append(t)
    return kept

def _category_plots_from_pheno_df(category_name: str,
                                  pheno_df: pd.DataFrame,
                                  selected_rows: dict[str, pd.Series],
                                  active_trait: str = None,         # ✅ 현재 선택된 trait
    active_category: str = None       # ✅ 현재 선택된 category
    ):
    """
    카테고리별 trait 플롯들을 grid 레이아웃으로 생성.
    - numeric: 전체 히스토그램 + 선택 품종 vertical line
    - categorical: 전체 bar count + 선택 품종 marker
    - date: interval binning + bar count + 선택 품종 marker
    """
    import plotly.graph_objects as go
    import numpy as np
    from dash import html, dcc
    traits_all = get_category_traits(category_name)
    valid_traits = _available_traits(pheno_df, traits_all, selected_rows)

    if not valid_traits:
        return html.Div([
            html.P(f"No valid traits found for {category_name}",
                   style={'textAlign': 'center', 'color': 'gray'})
        ], style={'width': '100%', 'padding': '20px'})

    colors = ['red', 'blue', 'green', 'purple', 'orange',
              'brown', 'pink', 'gray', 'olive', 'cyan']

    # --- date helper ---
    def create_date_intervals_optimized(values):
        if len(values) == 0:
            return [], [], []
        def get_interval_key(value):
            try:
                val = int(value)
                month, day = val // 100, val % 100
                if not (1 <= month <= 12 and 1 <= day <= 31):
                    return None
                if day <= 10: return f"{month}/1-10"
                elif day <= 20: return f"{month}/11-20"
                else: return f"{month}/21-last"
            except: return None
        def get_interval_range(key):
            month, day_range = key.split('/')
            month = int(month)
            if day_range == "1-10": return (float(f"{month}01"), float(f"{month}10"))
            elif day_range == "11-20": return (float(f"{month}11"), float(f"{month}20"))
            else: return (float(f"{month}21"), float(f"{month}31"))
        existing = {get_interval_key(v) for v in values if get_interval_key(v)}
        sorted_keys = sorted(existing, key=lambda x: (
            int(x.split('/')[0]),
            0 if x.endswith('1-10') else 1 if x.endswith('11-20') else 2
        ))
        ranges, labels, tips = [], [], []
        for k in sorted_keys:
            r = get_interval_range(k)
            if r:
                ranges.append(r)
                labels.append(k)
                tips.append(f"Month {k}")
        return ranges, labels, tips

    def find_value_interval(value, date_ranges):
        for i, (start, end) in enumerate(date_ranges):
            if start <= value <= end:
                return i
        return -1

    # --- grid layout ---
    plot_rows = []
    for i in range(0, len(valid_traits), 3):
        row_traits = valid_traits[i:i+3]
        row_plots = []

        for trait in row_traits:
            fig = go.Figure()
            series = pheno_df[trait]
            clean_data = series[series.notna()]

            selected_values = [(vid, row[trait]) for vid, row in selected_rows.items()
                               if trait in row and pd.notna(row[trait])]

            # ✅ fallback 함수로 trait type 판정
            trait_type = determine_trait_type_fallback(trait)

            # --- numeric ---
            if trait_type == 'numeric':
                vals = pd.to_numeric(clean_data, errors='coerce').dropna()
                if not vals.empty:
                    fig.add_trace(go.Histogram(
                        x=vals, nbinsx=30,
                        name="All Varieties",
                        marker_color="lightblue",
                        opacity=0.7,
                        showlegend=False
                    ))
                    hist_counts, _ = np.histogram(vals, bins=30)
                    max_count = max(hist_counts) if len(hist_counts) else 1
                else:
                    max_count = 1

                for j, (vid, val) in enumerate(selected_values):
                    try:
                        v = float(val)
                        c = colors[j % len(colors)]
                        fig.add_vline(
                            x=v,
                            line_dash='dash',
                            line_color=c,
                            line_width=2,
                            annotation_text=vid,
                            annotation_position='top',
                            annotation_font_size=10,
                            annotation_font_color=c
                        )
                    except: 
                        continue

            # --- categorical ---
            elif trait_type == 'categorical':
                vc = clean_data.astype(str).value_counts()
                fig.add_trace(go.Bar(
                    x=vc.index.astype(str),
                    y=vc.values,
                    name="All Varieties",
                    marker_color="lightcoral",
                    opacity=0.7,
                    showlegend=False
                ))
                max_count = vc.max() if not vc.empty else 1
                for j, (vid, val) in enumerate(selected_values):
                    cat = str(val)
                    if cat in vc.index:
                        c = colors[j % len(colors)]
                        fig.add_trace(go.Scatter(
                            x=[cat],
                            y=[max_count*1.1],
                            mode="markers+text",
                            marker=dict(color=c, size=12, symbol="diamond",
                                        line=dict(color="white", width=2)),
                            text=[vid],
                            textposition="top center",
                            textfont=dict(size=9, color=c),
                            name=vid
                        ))

            # --- date ---
            elif trait_type == 'date':
                vals = pd.to_numeric(clean_data, errors="coerce").dropna()
                if not vals.empty or selected_values:
                    date_ranges, date_labels, tooltips = create_date_intervals_optimized(vals)
                    counts = [len(vals[(vals >= s) & (vals <= e)]) for (s, e) in date_ranges]
                    fig.add_bar(
                        x=date_labels, y=counts,
                        marker_color="rgba(34,139,34,0.6)",
                        opacity=0.7,
                        customdata=tooltips,
                        hovertemplate='<b>%{customdata}</b><br>Count: %{y}<extra></extra>'
                    )
                    max_count = max(counts) if counts else 1
                    for j, (vid, val) in enumerate(selected_values):
                        try:
                            v = float(val)
                            idx = find_value_interval(v, date_ranges)
                            if idx >= 0:
                                c = colors[j % len(colors)]
                                fig.add_trace(go.Scatter(
                                    x=[date_labels[idx]],
                                    y=[max_count*1.1],
                                    mode="markers+text",
                                    marker=dict(color=c, size=12, symbol="diamond",
                                                line=dict(color="white", width=2)),
                                    text=[vid],
                                    textposition="top center",
                                    textfont=dict(size=9, color=c),
                                    name=vid
                                ))
                        except: 
                            continue
                    fig.update_xaxes(tickangle=45, automargin=True)

            fig.update_layout(
                margin=dict(l=10, r=10, t=10, b=10),
                height=280,
                showlegend=False
            )

            graph_item = html.Div([
                html.Div([
                    html.Button(
                        trait.replace('_', ' ').title(),
                        id={'type': 'trait-title-button', 'trait': trait, 'category': category_name},
                        key=f"{trait}-{category_name}-{active_trait}-{active_category}",  # ✅ key 추가!
                        className=(
                            "btn btn-outline-primary trait-title-btn"
                            if active_trait != trait or active_category != category_name
                            else "btn btn-outline-info trait-title-btn"
                        ),
                        style={
                            'padding': '6px 8px',
                            'fontSize': '12px',
                            'fontWeight': 'bold',
                            #'border': '2px solid #007bff',
                            'borderRadius': '6px',
                            'marginBottom': '8px',
                            'cursor': 'pointer',
                            'display': 'inline-block',
                            'whiteSpace': 'nowrap'
                        },
                        n_clicks=0,
                        title=f"Click to apply {trait.replace('_', ' ').title()} color mapping to pedigree"
                    )
                ], style={'textAlign': 'center', 'marginBottom': '8px'}),

                dcc.Graph(
                    figure=fig,
                    className="trait-individual-plot",
                    style={'width': '100%', 'height': '280px'},
                    config={'displayModeBar': False, 'webGL': False}
                )
            ],
            className="trait-graph-container",
            style={
                'width': '32%',
                'display': 'inline-block',
                'verticalAlign': 'top',
                'marginRight': '1%',
                'marginBottom': '15px',
                'border': '1px solid #e0e0e0',
                'borderRadius': '8px',
                'padding': '10px',
                'boxSizing': 'border-box'
            })
            row_plots.append(graph_item)

        plot_rows.append(html.Div(row_plots, style={'width': '100%', 'marginBottom': '10px', 'textAlign': 'left'}))

    return html.Div(plot_rows, className="trait-plots-container", style={'width': '100%'})



def create_phenotype_panel_unified(selected_nodes: list[str] | None,
                                   pheno_df: pd.DataFrame,
                                   table_enabled: bool = True):
    """
    - selected_nodes: IT number 리스트(문자열). 계보/무계보 모두 동일 규격.
    - pheno_df: PHENO_DF (전역에서 로드된 DataFrame)
    - table_enabled: 상단 DataTable 표시 여부
    """
   
    selected_nodes = selected_nodes or []
    selected_rows = _pick_rows_by_ids(pheno_df, selected_nodes)  # dict[id] = row(Series)
    #print(f"selected_rows: {selected_rows}")
    all_category_mappings = get_all_category_mappings()

    # 카테고리 탭 만들기
    tabs = []
    for category_name in all_category_mappings.keys():
        tab_content = html.Div([
            dcc.Loading(
                id=f'phenotype-loading-{category_name}',
                type='default',
                children=[
                    html.Div(
                        id=f'trait-spans-{category_name}',
                        className="mb-3"
                    ),
                    html.Div(
                        id=f'phenotype-plot-{category_name}',
                        children=_category_plots_from_pheno_df(category_name, pheno_df,  selected_rows or {},None,None),
                        style={'width': '100%'}
                    )
                ]
            )
        ], className="p-3")

        tabs.append(
            dbc.Tab(
                label=category_name.replace('_', ' ').title(),
                tab_id=f"phenotype-{category_name}",
                children=tab_content
            )
        )

    # 상단 테이블(옵션)


    # 범례(초기 숨김)
    '''
    legend_block = html.Div([
        html.H6("Trait Color Mapping", className="mb-2"),
        html.Div(id="trait-legend-container", className="mb-3")
    ], id="trait-legend-section", style={
        'padding': '10px',
        'backgroundColor': '#f8f9fa',
        'borderRadius': '8px',
        'border': '1px solid #dee2e6',
        'marginBottom': '15px',
        'width': '80%',
        'display': 'none'
    })
'''
    '''dbc.Row([
            dbc.Col([
                dbc.Button("📥 Download Phenotype Data CSV", id="download-phenotype-csv", color="primary", className="mb-2"),
                dcc.Download(id="download-phenotype")
            ], width=12)
        ])'''

    return html.Div([
        html.H5("Phenotype Distribution Analysis", className="mb-3"),
        #*table_block,
        #legend_block,
        dbc.Tabs(
            id="phenotype-category-tabs",
            active_tab="phenotype-yield_quality",
            children=tabs
        ),
        html.Hr(),
        
    ], className="p-3")

def _df_to_records(df: pd.DataFrame):
    return [] if df is None or df.empty else df.to_dict(orient='records')

def _records_to_df(records):
    return pd.DataFrame(records or [])

def _norm_combo_key(samples: List[str]) -> str:
    return "|".join(sorted({s for s in samples if s}))

def _to_jsonable_df(df: pd.DataFrame) -> List[Dict]:
    # dcc.Store 등에 넣을 때
    return df.to_dict(orient="records")

def _from_jsonable_df(obj) -> pd.DataFrame:
    return pd.DataFrame(obj or [])

def _build_variant_only(df: pd.DataFrame, sample_name: Optional[str]=None) -> pd.DataFrame:
    """
    Bar 차트용: Variant ID 유일 기준.
    필요한 컬럼: 'Variant ID','Chromosome','Position','Minor Allele','Trait','Subtrait' + 샘플 GT
    """
    use_cols = ['Variant ID','Chromosome','Position','Minor Allele','Trait','Subtrait']
    for c in use_cols:
        if c not in df.columns:
            df[c] = None

    out = df[use_cols + (['Allele'] if 'Allele' in df.columns else [])].copy()

    if sample_name and 'Allele' in out.columns:
        out.rename(columns={'Allele': f'{sample_name}_GT'}, inplace=True)

    out = out.sort_values(['Variant ID']).drop_duplicates(subset=['Variant ID'], keep='first')
    return out

def _build_variant_pvalue(df: pd.DataFrame, sample_name: Optional[str]=None) -> pd.DataFrame:
    """
    Scatter용: (Variant ID, P-value) 조합 유지.
    필요한 컬럼: 'Variant ID','Chromosome','Position','Minor Allele','MAF','P-value','PMID','Trait','Subtrait' + 샘플 GT
    """
    cols = ['Variant ID','Chromosome','Position','Minor Allele','MAF','P-value','PMID','Trait','Subtrait']
    for c in cols:
        if c not in df.columns:
            df[c] = None

    out = df[cols + (['Allele'] if 'Allele' in df.columns else [])].copy()

    if sample_name and 'Allele' in out.columns:
        out.rename(columns={'Allele': f'{sample_name}_GT'}, inplace=True)

    out = out.sort_values(['Variant ID','P-value']).drop_duplicates(subset=['Variant ID','P-value'], keep='first')
    return out

def _safe_error(requested, msg, combo_key=""):
    return {
        "status": "error",
        "requested_samples": list(requested or []),
        "matched_samples": [],
        "variant_only": [],
        "variant_pvalue": [],
        "error": str(msg),
        "combo_key": combo_key
    }

def _safe_empty(requested, msg="empty", combo_key=""):
    return {
        "status": "empty",
        "requested_samples": list(requested or []),
        "matched_samples": [],
        "variant_only": [],
        "variant_pvalue": [],
        "error": str(msg),
        "combo_key": combo_key
    }


# Local modules - 경로 문제 해결
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from pedi_module import PedigreeApp
    from gwas_module import GWASModule
    #from mannequin_callbacks import *
    #from mannequin_copy_callbacks import *
    GWAS_AVAILABLE = True
    print("✅ pedi_module과 gwas_module을 성공적으로 import했습니다.")
except ImportError as e:
    print(f"⚠️ GWAS module import failed: {e}")
    GWAS_AVAILABLE = False

# Initialize Cytoscape layouts
cyto.load_extra_layouts()

# VCF 정보 로드 함수


def load_phenotype_info():
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    
    
    pheno_path = os.path.join(BASE_DIR, "data","phenotype_data_en.csv")
    df=pd.read_csv(pheno_path, dtype=str)
    print(f"✅ Loaded {len(df)} phenotype data")
    print(f"✅ Phenotype data columns: {df.columns}")
    print(f"✅ Phenotype data columns: {df.head()}")
    return df

# 전역 로드
PHENO_DF = load_phenotype_info()





# Trait 정보 로드 함수들
def load_trait_info_data():
    """Excel 파일들로부터 trait 정보를 로드하고 concat"""
    excel_files = [
        '/data/vigor.xlsx', '/data/yield.xlsx', '/data/stress.xlsx', '/data/Plant_Quality.xlsx',
        '/data/plant_morphology.xlsx', '/data/groth_and_dev.xlsx', '/data/biological_process.xlsx', '/data/biochemical.xlsx'
    ]
    
    all_dfs = []
    
    for file in excel_files:
        try:
            file_path = os.path.join(os.path.dirname(__file__), file)
            df = pd.read_excel(file_path, skiprows=1)
            
            required_cols = ["Trait", "Trait ID", "Mapped Terms", "Description"]
            if all(col in df.columns for col in required_cols):
                df_selected = df[required_cols].copy()
                df_selected['Source'] = file.replace('.xlsx', '')
                df_selected['Group'] = file.replace('.xlsx', '')  # Add Group column
                all_dfs.append(df_selected)
                #print(f"✅ Loaded {len(df_selected)} traits from {file}")
                
        except Exception as e:
            print(f"⚠️ Error reading {file}: {e}")
    
    if all_dfs:
        trait_table = pd.concat(all_dfs, ignore_index=True)
        trait_table = trait_table.drop_duplicates(subset=['Trait ID'], keep='first')
        print(f"📊 Total unique traits loaded: {len(trait_table)}")
        return trait_table
    else:
        return pd.DataFrame()

def get_group_color_and_emoji(group_num: int):
    pastel_colors = px.colors.qualitative.Pastel1
    group_color = pastel_colors[group_num - 1] if group_num <= len(pastel_colors) else '#999'
    circle_symbol = html.Span("●", style={'color': group_color, 'font-size': '14px', 'marginRight': '6px'})
    return group_color, circle_symbol

def create_group_block(group_num, visible=True, removable=True, group_type="trait"):
    import plotly.express as px

    # 🎨 색상 리스트
    pastel_colors = px.colors.qualitative.Pastel1

    # ✅ 그룹 타입별 offset 계산
    if group_type == "trait":
        offset = 0
    elif group_type == "sample":
        offset = 5
    elif group_type == "variant":
        offset = 10
    else:
        offset = 0

    # ✅ 상대 인덱스 (1~5) → 동일한 색상 세트 공유
    relative_idx = group_num - offset
    color = pastel_colors[(relative_idx - 1) % len(pastel_colors)]

    # ✅ 라벨 및 dropdown 구성
    children = [
        html.Label(
            f"Group {relative_idx}:",
            style={
                'fontWeight': 'bold',
                'color': color,           # ✅ offset 반영 색상 적용
                'marginRight': '6px',
                'minWidth': '70px'
            }
        ),
        dcc.Dropdown(
            id={'type': 'group-dropdown', 'group_type': group_type, 'index': group_num},
            options=[],
            multi=True,
            placeholder="Select items...",
            className="group-dropdown",
            style={
                'flex': '1',
                'minWidth': '180px',
                'maxWidth': '240px',
                'fontSize': '13px'
            }
        )
    ]

    # ✖ 버튼 추가 (removable=True일 때만)
    if removable:
        children.append(
            html.Button(
                "✖",
                id={'type': 'remove-group-btn', 'group_type': group_type, 'index': group_num},
                style={
                    'border': 'none',
                    'background': 'none',
                    'color': color,       # ✅ 버튼 색도 동일 색상 사용
                    'cursor': 'pointer',
                    'fontSize': '14px'
                }
            )
        )

    # ✅ 그룹 block 스타일
    return html.Div(
        children,
        id={'type': 'group-block', 'group_type': group_type, 'index': group_num},
        style={
            'display': 'flex' if visible else 'none',
            'alignItems': 'center',
            'gap': '8px',
            'border': f'1px solid {color}',     # ✅ 그룹별 color 테두리
            'padding': '6px 8px',
            'borderRadius': '6px',
            'backgroundColor': '#f8f9fa',       # ✅ 카드 톤 통일
            'boxShadow': '0 1px 2px rgba(0,0,0,0.04)',
            'transition': 'all 0.3s ease',
            'flex': '0 1 calc(50% - 10px)',
            'minWidth': '300px',
        }
    )



# Trait 정보 전역 로드
TRAIT_INFO_DF = load_trait_info_data()

'''
def updated_check_sample_files_exist(entry_nos, data_path="/mnt/d/nabic/nabic2_6/fin_folder/filtered_results_2/matched_results_yesonly_2/"):
    """
    기존 check_sample_files_exist 함수의 개선 버전
    _SNP_gwas_yes.tsv 파일 형식에 특화
    """
    existing_samples = []
    
    if not os.path.exists(data_path):
        print(f"Data path not found: {data_path}")
        return existing_samples
    
    try:
        available_files = os.listdir(data_path)
        
        for entry_no in entry_nos:
            # entry_no_SNP_gwas_yes.tsv 형식의 파일 찾기
            target_file = f"{entry_no}_SNP_gwas_yes.tsv"
            if target_file in available_files:
                existing_samples.append(entry_no)
                
    except Exception as e:
        print(f"Error checking sample files: {e}")
        
    return existing_samples
'''

# =============================================================================
# DATABASE IMPORTS AND CONFIGURATION (pedi_utils.py 사용)
# =============================================================================

# DB 설정을 pedi_utils.py에서 가져오기
try:
    from assets.pedi_utils import create_header, create_footer, COMMON_STYLES  # 👈 올바른 경로
    DB_AVAILABLE = True
    
except ImportError as e:
    print(f"⚠️ pedi_utils.py import failed: {e}")
    # 백업: 직접 연결


# =============================================================================
# UNIFIED STYLES (ORDER2.MD requirement)
# =============================================================================

# =============================================================================
# HEADER AND FOOTER COMPONENTS (unified across app7_home.py, app7_1.py)
# =============================================================================

# =============================================================================
# URL PARAMETER PROCESSING (ORDER2.MD requirement)
# =============================================================================
def parse_url_params(search_string):
    """Parse URL parameters from search string"""
    if not search_string or not search_string.startswith('?'):
        return {}
    
    try:
        params = {}
        param_pairs = search_string[1:].split('&')
        for pair in param_pairs:
            if '=' in pair:
                key, value = pair.split('=', 1)
                params[key] = unquote(value)
        return params
    except Exception as e:
        print(f"Error parsing URL parameters: {e}")
        return {}

# =============================================================================
# GWAS 통합 기능 (test_apppedi.py에서 가져옴)
# =============================================================================
'''
def extract_vcf_values_with_fixed(app, selected_nodes, fixed_vcf_value):
    """VCF 값을 추출하되 고정값을 첫 번째로 우선배치"""
    if not selected_nodes:
        return [fixed_vcf_value] if fixed_vcf_value else [], {}
    
    # 기존 VCF 값들 추출
    all_vcf_values, vcf_info = extract_vcf_values_from_selected(app, selected_nodes)
    
    # 고정값을 첫 번째로 배치
    final_vcf_values = []
    if fixed_vcf_value:
        final_vcf_values.append(fixed_vcf_value)
        # 고정값을 제외한 나머지 값들 추가
        for value in all_vcf_values:
            if value != fixed_vcf_value and value not in final_vcf_values:
                final_vcf_values.append(value)
    else:
        final_vcf_values = all_vcf_values
    
    print(f"🔒 VCF 값 (고정값 우선): {final_vcf_values}")
    return final_vcf_values, vcf_info
'''

'''
def extract_vcf_values_from_selected(app, selected_nodes):
    """선택된 노드들에서 VCF 값들을 추출"""
    vcf_values = []
    vcf_info = []
    
    print(f"\n🔍 선택된 {len(selected_nodes)}개 품종에서 VCF 값 추출 중...")
    
    for node in selected_nodes:
        try:
            variety_info = app.get_variety_info_by_pedigree(str(node))
            # VCF Status와 status 필드 모두 확인 (호환성)
            vcf_status = variety_info.get("VCF Status") or variety_info.get("status")
            it_number = variety_info.get("IT Number")
            
            if vcf_status and vcf_status != 'No VCF':
                vcf_values.append(vcf_status)
                vcf_info.append({
                    'variety': node,
                    'it_number': it_number,
                    'vcf_value': vcf_status
                })
                print(f"   ✅ {node}: VCF = {vcf_status}, IT = {it_number}")
            else:
                print(f"   ❌ {node}: VCF 없음 (상태: {vcf_status})")
                
        except Exception as e:
            print(f"   ⚠️ {node}: VCF 확인 실패 - {e}")
            continue
    
    print(f"\n📋 추출된 VCF 값들: {vcf_values}")
    return vcf_values, vcf_info
'''
'''
def get_gwas_data_for_samples(vcf_values):
    """VCF 값들에 해당하는 GWAS 데이터를 가져옴 (서버 시작 안함)"""
    if not GWAS_AVAILABLE:
        print("❌ GWAS 모듈을 사용할 수 없습니다.")
        return None
        
    if not vcf_values:
        print("❌ VCF 값이 없어서 GWAS 데이터를 가져올 수 없습니다.")
        return None
    
    print(f"\n🧬 GWAS 데이터 로딩 중: {vcf_values}")
    
    try:
        # GWASModule 인스턴스 생성 (서버 시작 안함)
        gwas_module = GWASModule()
        
        # 전체 GWAS 데이터 가져오기
        gwas_data = gwas_module.get_gwas_data()
        available_samples = gwas_module.get_samples()
        
        print(f"   전체 GWAS 샘플: {len(available_samples)}개")
        print(f"   전체 GWAS 데이터: {len(gwas_data)}행")
        
        # 요청된 VCF 값들과 매칭되는 샘플들 찾기
        matched_samples = [s for s in vcf_values if s in available_samples]
        print(f"   매칭된 샘플: {matched_samples}")
        
        if not matched_samples:
            print(f"   ⚠️ 요청된 VCF 값들과 매칭되는 GWAS 샘플이 없습니다.")
            print(f"   요청된: {vcf_values}")
            print(f"   사용가능: {available_samples[:10]}...")  # 처음 10개만 표시
            return None
        
        # 매칭된 샘플들에 대한 데이터 필터링
        filtered_data = {}
        combined_sample_data = []
        
        for sample in matched_samples:
            sample_data = gwas_module.get_sample_data(sample)
            if not sample_data.empty:
                filtered_data[sample] = sample_data
                combined_sample_data.append(sample_data)
                print(f"   ✅ {sample}: {len(sample_data)}행의 GWAS 데이터")
        
        # 선택된 샘플들의 데이터를 합쳐서 full_data로 사용
        if combined_sample_data:
            import pandas as pd
            selected_samples_data = pd.concat(combined_sample_data, ignore_index=True)
            print(f"   📊 선택된 샘플들의 합쳐진 데이터: {len(selected_samples_data)}행")
            print(f"   📊 선택된 샘플들의 합쳐진 데이터컬럼: {selected_samples_data.columns}")
        else:
            selected_samples_data = pd.DataFrame()
        
        return {
            'gwas_module': gwas_module,
            'sample_data': filtered_data,  # 샘플별 개별 데이터
            'matched_samples': matched_samples,
            'full_data': selected_samples_data  # 선택된 샘플들의 합쳐진 데이터
        }
        
    except Exception as e:
        print(f"❌ GWAS 데이터 로딩 중 오류: {e}")
        return None
'''

def get_pedigree_variety_options():
    """
    계보도 CSV(=PedigreeApp)에서 사용 가능한 품종명 목록 조회 (자동완성용)
    DB 연결/쿼리 없음!
    """
    try:
        pedigree_app = get_pedigree_app()
        if not pedigree_app:
            return []

        # PedigreeApp 내부 CSV에서 cultivar_name 목록을 가져오는 메서드 사용
        # (이 메서드는 여러분의 PedigreeApp에 이미 구현되어 있어야 합니다:
        #   list_all_cultivars() -> List[str]
        # )
        names = pedigree_app.list_all_cultivars()
        names = [n for n in names if n and isinstance(n, str)]
        names = sorted(set(names))

        options = [{'label': n, 'value': n} for n in names]
        print(f"🌾 계보도 품종 옵션 로드(From PedigreeApp/CSV): {len(options)}개")
        return options

    except Exception as e:
        print(f"❌ get_pedigree_variety_options 오류: {e}")
        return []


def create_gwas_analysis_content(fixed_vcf_value):
    """GWAS 분석 탭 내용 생성 (고정값 포함)"""
    try:
        print(f"🔍 GWAS 분석 콘텐츠 생성: fixed_vcf_value={fixed_vcf_value}")
        
        if not fixed_vcf_value or fixed_vcf_value == '-':
            return html.Div([
                dbc.Alert([
                    html.H5("❌ No VCF data", className="alert-heading"),
                    html.P("VCF data is required for GWAS analysis.")
                ], color="danger")
            ])
        
        # 고정값 표시
        fixed_notice = f" (Fixed sample: '{fixed_vcf_value}')"
        
        return html.Div([
            
            
        html.Div(
            id='gwas-variant-section',
            style={
                'display': 'block',
                'margin-bottom': '25px',
                'borderBottom': '1px solid #e0e0e0',
                'paddingBottom': '8px'
            },
            children=[
                  dbc.Row(
                        [
                            # 왼쪽: 제목 + 설명
                            dbc.Col(
                            html.Div([
                                html.Div([
                                    html.H6(
                                        "Filter Apply & Trait Barplot",
                                        style={
                                            'marginBottom': '0px',
                                            'fontWeight': 'bold',
                                            'color': '#2c3e50'
                                        }
                                    ),
                                    html.I(
                                        "❓",
                                        id="trait-barplot-help-icon",
                                        style={
                                            "marginLeft": "1px",
                                            "cursor": "pointer",
                                            "color": "#6c757d",
                                            "fontSize": "16px",
                                            "position": "relative",
                                            "top": "-1px",      # ✅ 미세 정렬 (가장 안정적)
                                        }
                                    ),
                                    dbc.Tooltip(
                                        [
                                            "Apply filters (Subtrait, P-value, MAF, SNP Presence, Group)",
                                            html.Br(),
                                            "to generate a trait barplot showing only traits that meet the selected conditions."
                                        ],
                                        target="trait-barplot-help-icon",
                                        placement="right",
                                        style={"fontSize": "0.9em", "maxWidth": "480px"}
                                    ),
                                ],
                                style={
                                    "display": "flex",
                                    "alignItems": "center",
                                    "marginBottom": "4px"
                                }),
                            ]),
                            width=True,
                            style={
                                "display": "flex",
                                "flexDirection": "column",
                                "justifyContent": "center"
                            }
                        ),

                            # 오른쪽: 접기 버튼
                            dbc.Col(
                                dbc.Button(
                                    [html.I(className="fas fa-chevron-up", style={"marginRight": "6px"}), "Hide Filters"],
                                    id="filter-collapse-btn",
                                    color="secondary",
                                    size="sm",
                                    style={
                                        "padding": "2px 8px",
                                        "fontSize": "12px"
                                    }
                                ),
                                width="auto",
                                style={
                                    "display": "flex",
                                    "alignItems": "center",
                                    "justifyContent": "flex-end"
                                }
                            )
                        ],
                        justify="between",
                        align="center",
                        style={"marginBottom": "10px"}
                    ),
                        
                
                html.Hr(style={'marginTop': '5px', 'marginBottom': '5px'}),
                
                html.Div([
                        


                
                html.Div([
                # Sample Variant Positions Chart
                
                    # Radio button toggle for data view mode
                    

                    # Multi-Sample Analysis Modal (keep existing)
                    dbc.Modal([
                        dbc.ModalHeader(dbc.ModalTitle("Analyze with Additional Samples")),
                        dbc.ModalBody([
                            html.Div([
                                html.H4("Selected Samples", style={'color': '#2c3e50', 'marginBottom': '10px'}),
                                dash_table.DataTable(
                                    id='mannequin-combined-table',
                                    columns=[{'id':'Entry_No.','name':'Entry_No.'}],
                                    data=[],
                                    **base_table_props
                                )
                            ], id='mannequin-combined-display', style={'border':'1px solid #27ae60','padding':'12px',
                               'backgroundColor':'#f8fff8','marginBottom': '15px'}),
                            
                            # Tables for sample selection (simplified)
                            html.Div([
                                html.H4("Applied Varieties", style={'color': '#2c3e50', 'marginBottom': '10px'}),
                                html.Div(id='Mannequin_table',
                                    children=dash_table.DataTable(
                                        id='mannequin-table-a',
                                        columns=[{'id':'Entry_No.','name':'Entry_No.'}],
                                        data=[],
                                        **base_table_props
                                    ))
                            ], style={'border': '1px solid #bdc3c7', 'padding': '14px', 'marginBottom': '16px'})
                        ])
                    ], id="multi-sample-modal", size="xl", is_open=False),

                            # Modular Filter Options Section
                                                



                            html.Div([
                                # Filter Options Header with Toggle Button
                                html.Div([
                                    dbc.Button([
                                        "Filter Options ",
                                        html.I(id='filter-options-icon', className="fas fa-chevron-right")
                                    ], 
                                    id='filter-options-toggle',
                                    color="link",
                                    size="sm",
                                    style={
                                        'fontSize': '1rem', 
                                        'fontWeight': 'bold',
                                        'textDecoration': 'none',
                                        'padding': '2px 8px',
                                        'border': 'none'
                                    }),
                                    html.Span(id='filter-summary2', style={
                                        'marginLeft': '8px', 
                                        'fontSize': '12px', 
                                        'color': '#666'
                                    })
                                ], style={'marginBottom': '5px'}),
                                
                                # Collapsible Filter Container
                                dbc.Collapse([
                                    html.Div([
                                        html.Div([
                                        dcc.Checklist(
                                            id='filter-radio',
                                            options=[
                                                {'label': 'Trait (Subtrait)', 'value': 'trait'},
                                                {'label': 'P-value', 'value': 'pvalue'},
                                                {'label': 'MAF', 'value': 'maf'},
                                                {'label': 'SNP Presence', 'value': 'snp_presence'},
                                                {'label': 'Unique Mode', 'value': 'unique'}  # ✅ 추가됨    
                                                
                                            ],
                                            value=[],
                                            inline=True,
                                            inputStyle={"marginRight": "3px"},
                                            labelStyle={'margin-right': '10px', 'margin-left': '3px'},
                                            #style={'marginBottom': '10px','marginLeft':'15px'}
                                        ) 
                                        ], style={"marginBottom": "5px"}),
                                        # 필터들을 감싸는 Grid Wrapper
                                        html.Div([
                                            # Trait Filter Section
                                            html.Div([
                                                html.H6("Subtraits Selection", style={'marginBottom': '10px'}),
                                                html.Label("Choose Subtraits:", style={'margin-bottom': '3px', 'color': '#2c3e50'}),
                                                dcc.Dropdown(
                                                    id='integrated-subtrait-dropdown',
                                                    placeholder='Select one or more subtraits',
                                                    multi=True,
                                                    style={'margin-bottom': '10px', 'width': '70%'}
                                                ),
                                                html.P("Default: Shows all subtraits in sample-summary-barplot",
                                                    style={'fontSize': '0.85em', 'color': 'gray', 'fontStyle': 'italic'}),
                                            ], id='trait-section', style={
                                                'display': 'none',
                                                'padding': '10px',
                                                'border': '1px solid #ccc',
                                                'marginBottom': '10px',
                                                'borderRadius': '4px',
                                                'backgroundColor': '#f8f9fa'
                                            }),

                                            # P-value Filter Section
                                            html.Div([
                                                html.H6("P-value Filter", style={'marginBottom': '10px'}),
                                                html.Label("-log₁₀(P) Cutoff:", style={'margin-bottom': '3px', 'color': '#2c3e50'}),
                                                dcc.Input(
                                                    id="pvalue-cutoff",
                                                    type="number",
                                                    placeholder="5",
                                                    value=5,
                                                    min=0,
                                                    max=20,
                                                    step=0.1,
                                                    style={'width': '120px', 'margin-right': '10px'}
                                                ),
                                                html.P("(significance threshold)", style={'fontSize': '0.85em', 'color': 'gray', 'marginTop': '2px'}),
                                                html.P([
                                                    "By default, the scatter plot is rendered with a cutoff value of 5.",
                                                    html.Br(),
                                                    "This can be adjusted by the user."
                                                ], style={'fontSize': '0.85em', 'color': 'gray', 'fontStyle': 'italic', 'marginTop': '5px'}),
                                            ], id='pvalue-section', style={
                                                'display': 'none',
                                                'padding': '10px',
                                                'border': '1px solid #ccc',
                                                'marginBottom': '10px',
                                                'borderRadius': '4px',
                                                'backgroundColor': '#f8f9fa'
                                            }),

                                            # MAF Filter Section
                                            html.Div([
                                                html.H6("MAF Filter", style={'marginBottom': '10px'}),
                                                html.Label("MAF Threshold:", style={'margin-bottom': '3px', 'color': '#2c3e50'}),
                                                dcc.Input(
                                                    id="maf-cutoff",
                                                    type="number",
                                                    placeholder="0.05",
                                                    value=0.05,
                                                    min=0,
                                                    max=0.5,
                                                    step=0.01,
                                                    style={'width': '100px', 'margin-right': '10px'}
                                                ),
                                                html.P("(filter variants below this frequency)", style={'fontSize': '0.85em', 'color': 'gray', 'marginTop': '2px'}),
                                                dbc.Switch(id="maf-enabled", label="Enable MAF filtering", value=False, style={'marginTop': '5px'})
                                            ], id='maf-section', style={
                                                'display': 'none',
                                                'padding': '10px',
                                                'border': '1px solid #ccc',
                                                'marginBottom': '10px',
                                                'borderRadius': '4px',
                                                'backgroundColor': '#f8f9fa'
                                            }),
                                            # UNIQUE Filter Section
                                            html.Div([
                                                html.H6("Unique Variant Filter", style={'marginBottom': '10px'}),
                                                html.Label("Enable Unique Mode:", style={'margin-bottom': '3px', 'color': '#2c3e50'}),
                                                dbc.Switch(
                                                    id="unique-mode-enabled",
                                                    label="Show sample-unique variants only (requires 2–5 samples)",
                                                    value=False,
                                                    style={'marginTop': '5px'}
                                                ),
                                                html.P(
                                                    "This option highlights variants unique to each sample across the selected group.",
                                                    style={'fontSize': '0.85em', 'color': 'gray', 'marginTop': '5px', 'fontStyle': 'italic'}
                                                ),
                                                html.Div(id="unique-mode-warning",
                                                        style={'marginTop': '8px', 'fontSize': '12px', 'color': 'red'})
                                            ], id='unique-section', style={
                                                'display': 'none',
                                                'padding': '10px',
                                                'border': '1px solid #ccc',
                                                'marginBottom': '10px',
                                                'borderRadius': '4px',
                                                'backgroundColor': '#f8f9fa'
                                            }),

                                            # SNP Presence Section
                                            html.Div([
                                                html.H6("SNP Presence Threshold", style={'marginBottom': '10px'}),
                                                html.Label("Min Present Count:", style={'margin-bottom': '3px', 'color': '#2c3e50'}),
                                                dcc.Input(
                                                    id="snp-occurrence-min-count",
                                                    type="number",
                                                    placeholder="1",
                                                    value=1,
                                                    min=1,
                                                    step=1,
                                                    style={'width': '100px', 'margin-right': '10px'}
                                                ),
                                                html.P("(filter variants that appear in at least N samples)",
                                                    style={'fontSize': '0.85em', 'color': 'gray', 'marginTop': '2px'}),
                                                html.P(id='snp-presence-status',
                                                    style={'fontSize': '0.85em', 'color': '#888', 'fontStyle': 'italic', 'marginTop': '5px'}),
                                                dbc.Switch(id="snp-occurrence-enabled",
                                                        label="Enable SNP Presence Threshold", value=False, style={'marginTop': '5px'})
                                            ], id='snp-presence-section', style={
                                                'display': 'none',
                                                'padding': '10px',
                                                'border': '1px solid #ccc',
                                                'marginBottom': '10px',
                                                'borderRadius': '4px',
                                                'backgroundColor': '#f8f9fa'
                                            }),


                                        ],
                                        style={
                                            "display": "grid",
                                            "gridTemplateColumns": "repeat(auto-fit, minmax(320px, 1fr))",  # ✅ 항상 2열
                                            "gap": "15px"
                                        }),

                                          # 닫혀 있어도 항상 표시
                                    ],
                                    style={'width': '90%', 'margin': '0 auto 0 15px'})  # ✅ Collapse 박스 가운데 정렬
                            

                            
                                ], id='filter-options-collapse', is_open=False),
                            
                            
                            ]),
                            
                            html.Hr(style={'marginTop': '5px', 'marginBottom': '5px'}),


                            
                    
                    # Group Selection Section
                        html.Div([
                            # --- Group Options Header ---
                            html.Div([
                                dbc.Button(
                                    ["Group Selection ", html.I(id='group-options-icon', className="fas fa-chevron-right")],
                                    id='group-options-toggle',
                                    color="link",
                                    size="sm",
                                    style={
                                        'fontSize': '1rem',
                                        'fontWeight': 'bold',
                                        'textDecoration': 'none',
                                        'padding': '2px 8px',
                                        'border': 'none'
                                    }
                                ),
                                    html.Span(
                                        id='new-trait-alert',
                                        style={
                                            'marginLeft': '10px',
                                            'color': '#e67e22',
                                            'fontWeight': 'bold',
                                            'fontSize': '12px',
                                            'opacity': '0',
                                            'transition': 'opacity 0.5s ease'
                                        }
                                    ),
                                    dcc.Interval(id='new-trait-timer', interval=3000, n_intervals=0, disabled=True),

                                    # --- Variant ---
                                    html.Span(
                                        id='new-variant-alert',
                                        style={
                                            'marginLeft': '10px',
                                            'color': '#16a085',
                                            'fontWeight': 'bold',
                                            'fontSize': '12px',
                                            'opacity': '0',
                                            'transition': 'opacity 0.5s ease'
                                        }
                                    ),
                                    dcc.Interval(id='new-variant-timer', interval=3000, n_intervals=0, disabled=True),

                                    # --- Sample ---
                                    html.Span(
                                        id='new-sample-alert',
                                        style={
                                            'marginLeft': '10px',
                                            'color': '#2980b9',
                                            'fontWeight': 'bold',
                                            'fontSize': '12px',
                                            'opacity': '0',
                                            'transition': 'opacity 0.5s ease'
                                        }
                                    ),
                                    dcc.Interval(id='new-sample-timer', interval=3000, n_intervals=0, disabled=True),
                                #html.Span(id='group-summary', style={'marginLeft': '8px', 'fontSize': '12px', 'color': '#666'})
                            ], style={'marginBottom': '5px'}),

                            # --- Collapse Wrapper ---
                            dbc.Collapse([

                                # Radio for Group Type
                                html.Div([
                                    dcc.RadioItems(
                                        id='group-type-radio',
                                        options=[
                                            {'label': 'Trait', 'value': 'trait'},
                                            {'label': 'Variant ID', 'value': 'variant'},
                                            {'label': 'Sample', 'value': 'sample'}
                                        ],
                                        value='trait',
                                        inline=True,
                                        inputStyle={'marginRight': '3px'},
                                        labelStyle={'margin-right': '10px', 'margin-left': '3px'},
                                        style={'marginBottom': '5px', 'marginLeft': '15px'}
                                    ),
                                ], style={'marginBottom': '10px'}),

                                html.Hr(style={'marginTop': '5px', 'marginBottom': '10px'}),

                                # --- Group Containers ---
                                html.Div([
                                    html.Div(id='trait-group-container', style={
                                        'display': 'block',
                                        'padding': '10px',
                                        'border': '2px solid #ccc',
                                        'borderRadius': '6px',
                                        'marginBottom': '10px',
                                        'backgroundColor': '#fff',
                                        'overflow': 'hidden',
                                        'transition': 'max-height 0.4s ease',
                                        'maxHeight': '500px'
                                    }),
                                    html.Div(id='variant-group-container', style={
                                        'display': 'none',
                                        'padding': '10px',
                                        'border': '2px solid #ccc',
                                        'borderRadius': '6px',
                                        'marginBottom': '10px',
                                        'backgroundColor': '#fff',
                                        'overflow': 'hidden',
                                        'transition': 'max-height 0.4s ease',
                                        'maxHeight': '500px'
                                    }),
                                    html.Div(id='sample-group-container', style={
                                        'display': 'none',
                                        'padding': '10px',
                                        'border': '2px solid #ccc',
                                        'borderRadius': '6px',
                                        'marginBottom': '10px',
                                        'backgroundColor': '#fff',
                                        'overflow': 'hidden',
                                        'transition': 'max-height 0.4s ease',
                                        'maxHeight': '500px'
                                    }),
                                
                                    # --- Group Switch (아래로 이동) ---
                                
                                ], style={'width': '90%', 'margin': '0 auto 0 15px'}),

                                
                            ], id='group-options-collapse', is_open=False,
                            style={'width': '100%', 'margin': '0 auto 0 15px'})
                        ]),
                        html.Hr(style={'marginTop': '5px', 'marginBottom': '5px'}),
                    # Sample View Section
                    
                    #html.Hr(style={'marginTop': '5px', 'marginBottom': '5px'}),
                    # Plot mode only (table option removed per requirements)
                    dbc.Button(
                        [
                            html.I(className="fas fa-times me-2"),  # ❌ 아이콘 (me-2: 오른쪽 여백)
                            "Leave Group"                           # ✅ 텍스트
                        ],
                        id="close-group-barplot-btn",
                        title="Close group barplot",
                        color="link",  # ✅ 배경 투명 + 텍스트형 버튼 (Bootstrap "link" 스타일)
                        style={
                            "top": "8px",
                            "right": "12px",
                            "border": "none",
                            "background": "transparent",
                            "color": "#666",
                            "fontSize": "14px",
                            "cursor": "pointer",
                            "display": "none",  # ⛔️ 기본 숨김
                        },
                    ),
                    # Subtrait-based plot container (replaces integrated-trait-bar-container)
                    html.Div([
                        # 📊 Filter 기반 Barplot
                        html.Div([
                            dcc.Loading(
                                id='barplot-loading-filter',
                                type='default',
                                children=[
                                    dcc.Graph(
                                        id='integrated-trait-barplot',
                                        config={
                                            "displayModeBar": True,
                                            "displaylogo": False,
                                            "modeBarButtonsToRemove": ["pan2d", "lasso2d", "select2d"],
                                        },
                                        style={"width": "100%", "height": "100%"}
                                    )
                                ]
                            )
                        ], id='subtrait-container-wrapper-filter', style={'display': 'block'}),

                        
                        # 📊 Group 기반 Barplot
                        html.Div([
                            # ✅ 닫기 버튼 (우측 상단 고정)
                            

                            dcc.Loading(
                                id='barplot-loading-group',
                                type='default',
                                children=[
                                    dcc.Graph(
                                        id='integrated-trait-barplot-group',
                                        config={
                                            "displayModeBar": True,
                                            "displaylogo": False,
                                            "modeBarButtonsToRemove": ["pan2d", "lasso2d", "select2d"],
                                        },
                                        style={"width": "100%", "height": "100%"}
                                    )
                                ]
                            )
                        ], id='subtrait-container-wrapper-group', 
                        style={"position": "relative", 'display': 'none'})
                    ])
                    # Hidden components to maintain IDs for callbacks  
                    
                    # Scatter plot section (keep existing structure)
                


                            ])


                            ,
                                    
               
#---
                      
                        # Sample Summary Collapse 섹션들
                      
                    ], id='integrated-filter-section', style={ 'margin-bottom': '10px', 'max-width': '100%'}),
                    
                    
                        html.Div([
                            html.Div([
                            html.Hr(),

                            # 상단 행: 제목 + ? + 배지
                            html.Div([
                                # 왼쪽: 제목 + 도움말 아이콘
                                html.Div([
                                    html.H6("Result GWAS related SNP", style={"margin": 0}),
                                    html.I(
                                        "❓",
                                        id="gwas-help-icon",
                                        style={
                                            "marginLeft": "1px",
                                            "cursor": "pointer",
                                            "color": "#6c757d",
                                            "fontSize": "16px",
                                            "position": "relative",
                                            "top": "-1px",      # ✅ 미세 정렬 (가장 안정적)
                                        }
                                    ),
                                    dbc.Tooltip(
                                        [
                                            "Scatterplot: separated by chr1~chr12, displayed by Subtrait, and shows sample-specific positions in unique mode.",
                                            html.Br(),
                                            "Table: supports download and filtering.",
                                            html.Br(),
                                            "*Group is independent of sample option.*"
                                        ],
                                        target="gwas-help-icon",
                                        placement="right",
                                        style={"fontSize": "0.9em", "maxWidth": "480px"}
                                    ),
                                ], style={"display": "flex", "alignItems": "center"}),

                                # 오른쪽: Selected Variants 배지
                                dbc.Badge(
                                    "Selected Variants (0)",
                                    id='gwas-selected-badge',
                                    color='primary',
                                    pill=True,
                                    style={"cursor": "pointer", "marginBottom": "0px"}
                                ),
                            ],
                            style={
                                "display": "flex",
                                "justifyContent": "space-between",
                                "alignItems": "center",
                                "marginBottom": "6px"
                            }),

                            # 아래쪽 내용
                            html.Div(id='gwas-mini-table-container', style={'display': 'none'}),
                            html.Hr(style={'marginTop': '5px', 'marginBottom': '5px'}),
                        ]),

                            # 🔹 Filter-based Scatter/Table (Tabs)
                            dcc.Loading(
                                id='gwas-scatter-loading',
                                type='default',
                                children=[
                                    html.Div(
                                        id='gwas-sample-scatter-container',
                                        children=[
                                            dcc.Tabs(
                                                id='gwas-filter-tabs',
                                                value='scatter',
                                                children=[
                                                    dcc.Tab(label='Filter-based Scatter', value='scatter'),
                                                    dcc.Tab(label='Filter-based Table', value='table')
                                                ],
                                                style={'marginLeft': '15px', 'marginBottom': '8px'}
                                            ),

                                            # --- Scatter plot wrapper ---
                                            html.Div(
                                                id='gwas-sample-scatter-wrapper',
                                                children=[
                                                    dcc.Graph(
                                                        id='gwas-sample-scatter',
                                                        config={
                                                            'toImageButtonOptions': {
                                                                'format': 'png',
                                                                'filename': 'gwas_scatter_plot',
                                                                'height': 600,
                                                                'width': 1000,
                                                                'scale': 1
                                                            },
                                                            'displayModeBar': True,
                                                            'displaylogo': False,
                                                            'webGL': True,
                                                            'plotlyServerURL': False
                                                        }
                                                    )
                                                ],
                                                style={'display': 'block'}
                                            ),

                                            # --- Table wrapper ---
                                            html.Div(
                                                id='gwas-sample-table-wrapper',
                                                children=[
                                                    html.Div(dcc.Download(id="gwas-table-download")),
                                                    html.Div(id="gwas-table-download-status",
                                                            style={"color": "red", "fontSize": "0.9em", "marginTop": "5px"}),
                                                    dbc.Button(
                                                        "Download CSV",
                                                        id="gwas-table-download-btn",
                                                        n_clicks=0,
                                                        color="primary",
                                                        outline=True,
                                                        className="mb-2",
                                                        style={"marginTop": "10px", "marginLeft": "10px"}
                                                    ),
                                                    dash.dash_table.DataTable(
                                                        id='gwas-sample-table',
                                                        columns=[],
                                                        data=[],
                                                        page_size=10,
                                                        sort_action="native",
                                                        filter_action="native",
                                                        #selected_rows=[],
                                                        #row_selectable="multi",
                                                        style_table={'overflowX': 'auto'},
                                                        style_cell={'textAlign': 'left', 'padding': '5px'},
                                                        style_header={'fontWeight': 'bold'},
                                                    ),
                                                    html.Div(id="gwas-table-status",
                                                            style={"color": "red", "fontSize": "0.9em", "marginTop": "5px"})
                                                ],
                                                style={'display': 'none'}
                                            )
                                        ],
                                        style={'display': 'none'}
                                    ),

                                    # 🔹 Group-based Scatter/Table (Tabs)
                                    html.Div(
                                        id='gwas-sample-scatter-group-container',
                                        children=[
                                            dcc.Tabs(
                                                id='gwas-group-tabs',
                                                value='Gscatter',
                                                children=[
                                                    dcc.Tab(label='Group-based Scatter', value='Gscatter'),
                                                    dcc.Tab(label='Group-based Table', value='Gtable')
                                                ],
                                                style={'marginLeft': '15px', 'marginBottom': '8px'}
                                            ),

                                            # 그룹 정보
                                            html.Div(
                                                id="group-info-display",
                                                style={
                                                    "fontSize": "12px",
                                                    "color": "#444",
                                                    "margin": "10px auto",
                                                    "padding": "6px",
                                                    "border": "1px solid #ddd",
                                                    "borderRadius": "4px",
                                                    "width": "90%",
                                                    "display": "none"
                                                }
                                            ),

                                            # --- Group Scatter ---
                                            html.Div(
                                                id='gwas-sample-scatter-group-wrapper',
                                                children=[
                                                    dcc.Graph(
                                                        id='gwas-sample-scatter-group',
                                                        config={
                                                            'displayModeBar': True,
                                                            'displaylogo': False,
                                                            'toImageButtonOptions': {
                                                                'format': 'png',
                                                                'filename': 'gwas_group_scatter',
                                                                'height': 600,
                                                                'width': 1000,
                                                                'scale': 1
                                                            },
                                                            'webGL': True
                                                        }
                                                    )
                                                ],
                                                style={'display': 'block'}
                                            ),

                                            # --- Group Table ---
                                            html.Div(
                                                id='gwas-sample-table-group-wrapper',
                                                children=[
                                                    html.Div(dcc.Download(id="gwas-table-download-group")),
                                                    dbc.Button(
                                                        "Download Group CSV",
                                                        id="gwas-table-download-btn-group",
                                                        n_clicks=0,
                                                        color="primary",
                                                        outline=True,
                                                        className="mb-2",
                                                        style={"marginTop": "10px", "marginLeft": "10px"}
                                                    ),
                                                    dash.dash_table.DataTable(
                                                        id='gwas-sample-table-group',
                                                        columns=[],
                                                        data=[],
                                                        page_size=10,
                                                        sort_action="native",
                                                        filter_action="native",
                                                        #selected_rows=[],
                                                        #row_selectable="multi",
                                                        style_table={'overflowX': 'auto'},
                                                        style_cell={'textAlign': 'left', 'padding': '5px'},
                                                        style_header={'fontWeight': 'bold'},
                                                    )
                                                ],
                                                style={'display': 'none'}
                                            )
                                        ],
                                        style={'display': 'none'}
                                    )
                                ],
                                style={'width': '100%'}
                            )
                        ], style={'margin-bottom': '20px'})
                



            ]),
            
            # 숨겨진 Store들 (상태 관리용)
           # dcc.Store(id='gwas-selected-samples-store', data=[fixed_vcf_value]),
          # dcc.Store(id='snp-occurrence-min-count-store', data=1),
           dcc.Store(id='mannequin-selected-entry-nos', data=[]), 
           dcc.Store(id='mannequin-selected-entry-nos-copy2', data=[]), 
                 # VCF_INFO_DF 데이터를 저장할 스토어
            dcc.Store(id='mannequin-active-filter', data=None), 
            dcc.Store(id='mannequin-active-filter-copy2', data=None), 
        #)
    ], className="p-3")
        
    except Exception as e:
        print(f"❌ Error creating GWAS content: {e}")
        return html.Div([
            dbc.Alert([
                html.H5("❌ Error", className="alert-heading"),
                html.P(f"Error: {str(e)}")
            ], color="danger")
        ])
# =============================================================================
# VARIETY INFO FUNCTIONS (누락된 핵심 함수들 추가)
# =============================================================================


@callback(
    Output('vcf-info-global', 'data'),
    Input('app-init', 'data'),  # 앱 초기화 시 실행
    prevent_initial_call=False
)
def load_vcf_data_to_store(_):
    """
    VCF_INFO_DF 데이터를 글로벌 store에 로드
    """
    global VCF_INFO_DF
    if not VCF_INFO_DF.empty:
        return VCF_INFO_DF.to_dict('records')
    return []


def get_variety_info(variety_id) -> Optional[Dict]:
    if not variety_id:
        return None

    vid = _norm_id(variety_id)
    try:
        rt = load_resource_types()
        vm = load_vcf_mapping()
        if 'id' in rt.columns:
            rt = rt.copy(); rt['id'] = rt['id'].map(_norm_id)
        if 'variety_id' in vm.columns:
            vm = vm.copy(); vm['variety_id'] = vm['variety_id'].map(_norm_id)

        row = rt.loc[rt['id'] == vid].head(1)

        # variety_type: phenotype.variety_type -> resource_type -> 'Unknown'
        variety_type = 'Unknown'
        if _HAS_PHENO:
            try:
                ph = load_phenotype_info()
                if 'id' in ph.columns:
                    ph = ph.copy(); ph['id'] = ph['id'].map(_norm_id)
                    hit = ph.loc[ph['id'] == vid]
                    if not hit.empty:
                        if 'variety_type' in hit.columns and pd.notna(hit.iloc[0]['variety_type']):
                            variety_type = str(hit.iloc[0]['variety_type']).strip() or 'Unknown'
                        elif 'resource_type' in hit.columns and pd.notna(hit.iloc[0]['resource_type']):
                            variety_type = str(hit.iloc[0]['resource_type']).strip() or 'Unknown'
            except Exception:
                pass
        if variety_type == 'Unknown' and not row.empty and 'resource_type' in row.columns:
            if pd.notna(row.iloc[0]['resource_type']):
                variety_type = str(row.iloc[0]['resource_type']).strip() or 'Unknown'

        vcf_map = vm.set_index('variety_id')['vcf_tag'] if not vm.empty else pd.Series(dtype=str)
        status = vcf_map.get(vid, None)
        status = status if (isinstance(status, str) and status.strip()) else 'No VCF'

        if row.empty:
            return {
                'id': vid,
                'variety_name_en': '-',
                'line_name_en': '-',
                'variety_type': variety_type if variety_type else 'Unknown',
                'status': status
            }

        r0 = row.iloc[0]
        variety_name_en = _dash_if_empty(r0['variety_name_en'] if 'variety_name_en' in row.columns else None)
        line_name_en    = _dash_if_empty(r0['line_name_en'] if 'line_name_en' in row.columns else None)

        return {
            'id': vid,
            'variety_name_en': variety_name_en,
            'line_name_en': line_name_en,
            'variety_type': variety_type if variety_type else 'Unknown',
            'status': status
        }

    except Exception as e:
        print(f"[get_variety_info CSV] error for {vid}: {e}")
        return {
            'id': vid,
            'variety_name_en': '-',
            'line_name_en': '-',
            'variety_type': 'Unknown',
            'status': 'No VCF'
        }

def create_header_section(variety_id):
    """품종 정보 카드 생성 (사용자 원본 함수)"""
    variety_data = get_variety_info(variety_id) if variety_id else None
    
    return html.Div([
        html.H2("Rice Variety Information Card", style={
            'textAlign': 'center',
            'color': '#2c3e50',
            'marginBottom': '20px',
            'fontSize': '2rem',
            'fontWeight': 'bold'
        }),
        html.Div([
            html.Table([
                html.Tr([
                    html.Td("Variety Number:", style={'fontWeight': 'bold', 'padding': '8px', 'backgroundColor': '#f8f9fa'}),
                    html.Td(variety_data['id'] if variety_data is not None else "-", style={'padding': '8px'}),
                    html.Td("Variety Name:", style={'fontWeight': 'bold', 'padding': '8px', 'backgroundColor': '#f8f9fa'}),
                    html.Td(variety_data['variety_name_en'] if variety_data is not None else "-", style={'padding': '8px'})
                ]),
                html.Tr([
                    html.Td("Line Name:", style={'fontWeight': 'bold', 'padding': '8px', 'backgroundColor': '#f8f9fa'}),
                    html.Td(variety_data['line_name_en'] if variety_data is not None else "-", style={'padding': '8px'}),
                    html.Td("Type:", style={'fontWeight': 'bold', 'padding': '8px', 'backgroundColor': '#f8f9fa'}),
                    html.Td(variety_data['variety_type'] if variety_data is not None else "-", style={'padding': '8px'})
                ]),
                html.Tr([
                    html.Td("Status:", style={'fontWeight': 'bold', 'padding': '8px', 'backgroundColor': '#f8f9fa'}),
                    html.Td(variety_data['status'] if variety_data is not None else "-", style={'padding': '8px'}),
                    html.Td("", style={'padding': '8px'}),
                    html.Td("", style={'padding': '8px'})
                ])
            ], style={
                'width': '100%',
                'borderCollapse': 'collapse',
                'border': '1px solid #dee2e6',
                'borderRadius': '8px',
                'overflow': 'hidden',
                'backgroundColor': 'white',
                'boxShadow': '0 2px 4px rgba(0,0,0,0.1)',
                'margin': 'auto',
                'maxWidth': '800px'
            })
        ], style={'margin': '20px'})
    ], style={
        'backgroundColor': '#ffffff',
        'padding': '20px',
        'marginBottom': '20px',
        'borderRadius': '8px',
        'boxShadow': '0 2px 4px rgba(0,0,0,0.1)'
    })

# =============================================================================
# GLOBAL INSTANCES (한 번만 생성)
# =============================================================================
PEDIGREE_APP_INSTANCE = None

# Pedigree 정보 캐시 (v2024 최적화)
_pedigree_info_cache = {}

# GWAS 데이터 캐시 (v2024 최적화)
_gwas_module_instance = None
_gwas_base_data = None  # 전체 GWAS 데이터 캐시
_gwas_available_samples = None  # 사용가능한 샘플 목록
_gwas_sample_cache = {}  # 개별 샘플 데이터 캐시

def get_pedigree_app():
    """전역 PedigreeApp 인스턴스 반환 (싱글톤 패턴)"""
    global PEDIGREE_APP_INSTANCE
    if PEDIGREE_APP_INSTANCE is None:
        PEDIGREE_APP_INSTANCE = PedigreeApp()
        print("✅ PedigreeApp 인스턴스 생성됨")
    return PEDIGREE_APP_INSTANCE

def get_variety_info_by_pedigree_cached(pedigree_name):
    """
get_variety_info_by_pedigree의 캐시된 버전 - 성능 최적화 (v2024)
동일한 pedigree_name에 대한 반복 호출을 방지하여 DB 부하 감소
"""
    global _pedigree_info_cache
    
    cache_key = str(pedigree_name)
    
    # 캐시 히트
    if cache_key in _pedigree_info_cache:
        print(f"📋 Pedigree cache hit for {cache_key}")
        return _pedigree_info_cache[cache_key].copy()
    
    # 캐시 미스 - 실제 데이터 조회
    pedigree_app = get_pedigree_app()
    try:
        variety_info = pedigree_app.get_variety_info_by_pedigree(pedigree_name)
        
        # 캐시에 저장 (최대 200개 제한)
        if len(_pedigree_info_cache) < 200:
            _pedigree_info_cache[cache_key] = variety_info.copy() if isinstance(variety_info, dict) else variety_info
            print(f"📦 Cached pedigree info for {cache_key}")
        
        return variety_info
    except Exception as e:
        print(f"❌ Error getting pedigree info for {cache_key}: {e}")
        return {}

def clear_pedigree_cache():
    """학습될 pedigree 캐시 클리어"""
    global _pedigree_info_cache
    _pedigree_info_cache.clear()
    print("🗑️ Pedigree cache cleared")


def get_gwas_module_singleton():
    global _gwas_module_instance
    if not GWAS_AVAILABLE:
        print("❌ GWAS not available")
        return None, None, None
    if _gwas_module_instance is None:
        try:
            print("🧬 Initializing GWAS module singleton...")
            _gwas_module_instance = GWASModule()  # ✅ 경로 명시
            samples = _gwas_module_instance.get_samples()
            print(f"✅ GWAS singleton ready: {len(samples)} samples")
        except Exception as e:
            print(f"❌ GWAS singleton initialization failed: {e}")
            return None, None, None
    return _gwas_module_instance, None, _gwas_module_instance.get_samples()

'''
def get_gwas_module_singleton():
    """
GWAS 모듈 싱글톤 인스턴스 - 기본 데이터 미리 로드 (v2024)
반복적인 GWASModule 생성을 방지하여 성능 향상
"""
    global _gwas_module_instance, _gwas_base_data, _gwas_available_samples
    
    if not GWAS_AVAILABLE:
        return None, None, None
    
    # 처음 호출 시 초기화
    if _gwas_module_instance is None:
        print("🧬 Initializing GWAS module singleton...")
        try:
            _gwas_module_instance = GWASModule()
            _gwas_base_data = _gwas_module_instance.get_gwas_data()
            _gwas_available_samples = _gwas_module_instance.get_samples()
            print(f"✅ GWAS singleton initialized: {len(_gwas_available_samples)} samples, {len(_gwas_base_data)} rows")
        except Exception as e:
            print(f"❌ GWAS singleton initialization failed: {e}")
            return None, None, None
    
    return _gwas_module_instance, _gwas_base_data, _gwas_available_samples
'''


'''
def get_gwas_data_for_samples_optimized(vcf_values, use_cache=True):
    """
최적화된 GWAS 데이터 로드 - checkbox-driven lazy loading (v2024)

파이: 
- 싱글톤 GWAS 모듈 인스턴스
- 개별 샘플 데이터 캐시
- 사용자 체크박스 선택 시만 데이터 로드
"""
    if not vcf_values:
        print("❌ No VCF values provided for GWAS")
        return None
    
    # GWAS 모듈 싱글톤 가져오기
    gwas_module, base_data, available_samples = get_gwas_module_singleton()
    if not gwas_module:
        print("❌ GWAS module not available")
        return None
    
    # 매칭되는 샘플 찾기
    matched_samples = [s for s in vcf_values if s in available_samples]
    if not matched_samples:
        print(f"⚠️ No matching GWAS samples for: {vcf_values}")
        return None
    
    print(f"🧬 Optimized GWAS loading: {len(matched_samples)} samples")
    
    # 캐시된 또는 새로 로드되는 샘플 데이터
    sample_data_dict = {}
    combined_sample_data = []
    cache_hits = 0
    
    for sample in matched_samples:
        if use_cache and sample in _gwas_sample_cache:
            # 캐시 히트
            sample_data = _gwas_sample_cache[sample]
            cache_hits += 1
        else:
            # 새로 로드
            sample_data = gwas_module.get_sample_data(sample)
            if use_cache and len(_gwas_sample_cache) < 50:  # 최대 50개 샘플 캐시
                _gwas_sample_cache[sample] = sample_data
        
        if not sample_data.empty:
            sample_data_dict[sample] = sample_data
            combined_sample_data.append(sample_data)
    
    print(f"📋 GWAS Cache performance: {cache_hits}/{len(matched_samples)} hits")
    
    # 결과 데이터 구성
    if combined_sample_data:
        import pandas as pd
        selected_samples_data = pd.concat(combined_sample_data, ignore_index=True)
        
        result = {
            'gwas_module': gwas_module,
            'full_data': selected_samples_data,
            'sample_data': sample_data_dict,
            'available_samples': available_samples,
            'requested_samples': vcf_values,
            'matched_samples': matched_samples
        }
        
        print(f"✅ GWAS data loaded efficiently: {len(selected_samples_data)} rows")
        return result
    
    return None
'''
# 전역 스토어(단 하나)
_GWAS_STORE = {"per_sample": {}, "combos": {}}

def _df_to_records(df):
    import pandas as pd
    try:
        return [] if df is None or getattr(df, 'empty', True) else df.to_dict(orient='records')
    except Exception as e:
        print(f"[df_to_records] ❌ {e}")
        return []

def _safe_error(requested, msg, combo_key=""):
    return {
        "status": "error",
        "requested_samples": list(requested or []),
        "matched_samples": [],
        "variant_only": [],
        "variant_pvalue": [],
        "error": str(msg),
        "combo_key": combo_key
    }

def _safe_empty(requested, msg="empty", combo_key=""):
    return {
        "status": "empty",
        "requested_samples": list(requested or []),
        "matched_samples": [],
        "variant_only": [],
        "variant_pvalue": [],
        "error": str(msg),
        "combo_key": combo_key
    }
def get_gwas_data_for_samples_optimized(
    vcf_values,
    use_cache=True,
    return_jsonable=True,
    external_store=None
):
    """
    선택된 샘플 ID 리스트(vcf_values)를 받아
    - per-sample 데이터(variant_only, variant_pvalue) 캐싱
    - 조합(combos) 데이터 구성

    variant_only: Barplot용 (Variant ID 기준 유일 + 샘플별 GT 열)
    variant_pvalue: Scatter용 (Variant ID, P-value 기준 + 샘플별 GT 열)
    """
    import time, pandas as pd

    try:
        print("\n[GWAS-LOADER] ▶ start")
        print(f"[GWAS-LOADER] vcf_values={vcf_values}")

        if not vcf_values:
            return _safe_empty(vcf_values, "no samples")

        store = external_store if isinstance(external_store, dict) else _GWAS_STORE

        gwas_module, _, available_samples = get_gwas_module_singleton()
        if not gwas_module:
            return _safe_error(vcf_values, "gwas module not available")

        matched = [s for s in vcf_values if s in (available_samples or [])]
        print(f"[GWAS-LOADER] matched={matched}")

        combo_key = "|".join(sorted(set(matched))) if matched else ""
        if not matched:
            return _safe_empty(vcf_values, "no matching samples", combo_key)

        combos = store.setdefault("combos", {})
        per_sample = store.setdefault("per_sample", {})

        # 캐시 히트
        if use_cache and combo_key in combos and combos[combo_key].get("status") == "ready":
            combo = combos[combo_key]
            return {
                "status": "ready",
                "requested_samples": combo["requested_samples"],
                "matched_samples": combo["matched_samples"],
                "variant_only": _df_to_records(combo["variant_only"]),
                "variant_pvalue": _df_to_records(combo["variant_pvalue"]),
                "error": combo.get("error", ""),
                "combo_key": combo_key,
            }

        # 로딩 상태 기록
        combos[combo_key] = {
            "status": "loading",
            "requested_samples": list(vcf_values),
            "matched_samples": list(matched),
            "variant_only": None,
            "variant_pvalue": None,
            "error": "",
            "ts": time.time()
        }

        # ── 1) per-sample 데이터 로드
        vo_list, vp_list = [], []
        for s in matched:
            try:
                sdf = gwas_module.get_sample_data(s)
                print(f"[GWAS-LOADER] sample={s} shape={getattr(sdf,'shape',None)}")

                vo = _build_variant_only(sdf, sample_name=s)
                vp = _build_variant_pvalue(sdf, sample_name=s)

                vo_list.append(vo)
                vp_list.append(vp)

                # per-sample 캐시에 저장
                per_sample[s] = {
                    "status": "ready",
                    "variant_only": vo,
                    "variant_pvalue": vp,
                    "error": ""
                }

            except Exception as e:
                combos[combo_key].update({"status": "error", "error": f"sample {s} load error: {e}"})
                return _safe_error(vcf_values, f"sample {s} load error: {e}", combo_key)

        # ── 2) 콤보 조합
        # Bar용
        if vo_list:
            vo_all = pd.concat(vo_list, ignore_index=True)
            base = vo_all[['Variant ID','Chromosome','Position','Minor Allele','Trait','Subtrait']].drop_duplicates(subset=['Variant ID'])
            variant_only_combo = base
            for vo in vo_list:
                gt_cols = [c for c in vo.columns if c.endswith("_GT")]
                variant_only_combo = variant_only_combo.merge(
                    vo[['Variant ID'] + gt_cols], on='Variant ID', how='left'
                )
        else:
            variant_only_combo = pd.DataFrame(columns=['Variant ID','Chromosome','Position','Trait','Subtrait'] + [f"{s}_GT" for s in matched])

        # Scatter용
        if vp_list:
            # dtype 강제 통일
            for i in range(len(vp_list)):
                vp_list[i]['Variant ID'] = vp_list[i]['Variant ID'].astype(str)
                if 'P-value' in vp_list[i].columns:
                    vp_list[i]['P-value'] = pd.to_numeric(vp_list[i]['P-value'], errors='coerce')
                if 'Position' in vp_list[i].columns:
                    vp_list[i]['Position'] = pd.to_numeric(vp_list[i]['Position'], errors='coerce')

            vp_all = pd.concat(vp_list, ignore_index=True)

            # base 만들기
            base = vp_all[['Variant ID','P-value','Chromosome','Position','Minor Allele','MAF','PMID','Trait','Subtrait']].drop_duplicates()

            variant_pvalue_combo = base
            for vp in vp_list:
                gt_cols = [c for c in vp.columns if c.endswith("_GT")]
                if gt_cols:
                    variant_pvalue_combo = variant_pvalue_combo.merge(
                        vp[['Variant ID','P-value'] + gt_cols].astype({
                            'Variant ID': str,
                            'P-value': float
                        }),
                        on=['Variant ID','P-value'],
                        how='left'
                    )
        else:
            variant_pvalue_combo = pd.DataFrame(columns=[
                'Variant ID','Chromosome','Position','Minor Allele','MAF','P-value','PMID','Trait','Subtrait'
            ])

        combos[combo_key].update({
            "variant_only": variant_only_combo,
            "variant_pvalue": variant_pvalue_combo,
            "status": "ready",
            "error": ""
        })

        print(f"[GWAS-LOADER] ✅ done combo_key={combo_key} "
              f"rows(vo/vp)=({len(variant_only_combo)}, {len(variant_pvalue_combo)})")

        return {
            "status": "ready",
            "requested_samples": combos[combo_key]["requested_samples"],
            "matched_samples": combos[combo_key]["matched_samples"],
            "variant_only": _df_to_records(variant_only_combo),
            "variant_pvalue": _df_to_records(variant_pvalue_combo),
            "error": "",
            "combo_key": combo_key,
        }

    except Exception as e:
        print(f"[GWAS-LOADER] ❌ fatal: {e}")
        return _safe_error(vcf_values, f"fatal: {e}")



def clear_gwas_cache():
    """학습될 GWAS 캐시 전체 클리어"""
    global _gwas_module_instance, _gwas_base_data, _gwas_available_samples, _gwas_sample_cache
    
    _gwas_module_instance = None
    _gwas_base_data = None
    _gwas_available_samples = None
    _gwas_sample_cache.clear()
    
    print("🗑️ GWAS cache completely cleared")

# =============================================================================
# MAIN LAYOUT FUNCTION (함수형으로 변경)
# =============================================================================
def create_layout(processed_name='', variety_id=''):
    """URL 접근 시 자동 GWAS 활성화"""
    
    # CSS 애니메이션 스타일 추가
    
    
    # 1. 품종 정보 카드 생성
    header_section = create_header_section(variety_id) if variety_id else html.Div()
    
    # 2. 검색값 설정
    search_value = ""
    
    # 3. 자동 pedigree 생성 + 분석 탭 초기화
    initial_pedigree_elements = []
    initial_analysis_tabs = [
        dbc.Tab(label="Phenotype Analysis", tab_id="phenotype-tab"),
        dbc.Tab(label="GWAS Analysis", tab_id="gwas-tab", disabled=True),
    ]
    initial_active_tab = "phenotype-tab"
    initial_vcf_values = []
    
    # VCF 상태 변수 초기화
    vcf_status = None
    fixed_vcf_item_data = None
    fixed_vcf_item_nopedi_data = None
    parsed_id = (
        variety_id if processed_name == 'nopedi'
        else processed_name
    )
    #print(processed_name)
    if processed_name and processed_name == 'nopedi':
        # VCF 전용 모드 - variety_id로 VCF 상태 확인 및 분석 탭 활성화
        if variety_id:
            variety_data = get_variety_info(variety_id)
            #print(variety_data)
            if variety_data:
                # variety_data의 'status' 키가 VCF 상태를 포함 (vcf_tag) - 안전한 처리
                if isinstance(variety_data, dict):
                    vcf_status = variety_data.get('status', 'No VCF')
                else:
                    vcf_status = 'No VCF'
                
                # nopedi case: 2 fields (variety_id, status) - 안전한 처리
                if isinstance(variety_data, dict):
                    fixed_vcf_item_nopedi_data = {
                        'variety_id': variety_data.get('id', variety_id),
                        'status': vcf_status
                    }
                    variety_name = variety_data.get('variety_name_en', '')
                else:
                    fixed_vcf_item_nopedi_data = {
                        'variety_id': variety_id,
                        'status': vcf_status
                    }
                    variety_name = str(variety_data) if variety_data else ''
                
                print(f"🔍 VCF status from variety_id: {vcf_status}")
                print(f"🔍 Variety name from variety_id: {variety_name}")

                
                # nopedi case: 항상 phenotype 분석 활성화, VCF가 있으면 GWAS도 활성화
                if variety_name:  # variety_name이 있으면 항상 phenotype 분석 가능
                    # VCF 여부에 따라 GWAS 탭 활성화/비활성화
                    if vcf_status and vcf_status != 'No VCF':
                        # VCF가 있는 경우: GWAS 탭도 활성화
                        initial_vcf_values = [vcf_status]
                        initial_analysis_tabs = [
                            dbc.Tab(label="Phenotype Analysis", tab_id="phenotype-tab"),
                            dbc.Tab(label=f"GWAS Analysis", tab_id="gwas-tab", disabled=False),
                        ]
                        print(f"✅ nopedi case: GWAS enabled with VCF: {vcf_status}")
                    else:
                        # VCF가 없는 경우: phenotype만 활성화, GWAS 탭 완전 제거
                        initial_analysis_tabs = [
                            dbc.Tab(label="Phenotype Analysis", tab_id="phenotype-tab"),
                        ]
                        print(f"✅ nopedi case: Only phenotype enabled (GWAS tab hidden - No VCF)")
                    
                    # VCF 여부와 관계없이 phenotype 분석 초기 콘텐츠 항상 생성
                    initial_active_tab = "phenotype-tab"
                    try:
                        if PHENO_DF is None or PHENO_DF.empty:
                            initial_tab_content = html.Div([
                                dbc.Alert([
                                    html.H5("❌ No Phenotype Data", className="alert-heading"),
                                    html.P("PHENO_DF is empty or not loaded.")
                                ], color="danger")
                            ])
                        else:
                            # variety_id(IT 번호)만 넘겨 통합 패널 생성
                            print(f"variety_id: {variety_id}")
                            initial_tab_content = create_phenotype_panel_unified(
                                selected_nodes=[str(variety_id)],
                                pheno_df=PHENO_DF,
                                table_enabled=True
                            )
                        print(f"✅ nopedi case: phenotype panel created for {variety_id}")
                    except Exception as e:
                        print(f"⚠️ nopedi phenotype panel creation failed: {e}")
                        initial_tab_content = html.Div([
                            dbc.Alert([
                                html.H5("⚠️ Content Creation Failed", className="alert-heading"),
                                html.P(f"Error: {str(e)}")
                            ], color="warning")
                        ])
        
    elif processed_name and processed_name != 'nopedi':
        # 🎯 자동 pedigree 생성 + GWAS 초기화
        search_value = processed_name
        try:
            pedigree_app = get_pedigree_app()

            nodes, edges = pedigree_app.get_connected_nodes(processed_name, 2, 2)

            initial_pedigree_elements = pedigree_app.create_cytoscape_elements(nodes, edges)
            for element in initial_pedigree_elements:
                if 'source' not in element['data'] and element['data']['id'] == str(processed_name):
                    element['classes'] = (element.get('classes','') + ' search-result').strip()
            
            print(f"initial_pedigree_elements: {initial_pedigree_elements}")
            variety_info = get_variety_info_by_pedigree_cached(processed_name) or {}
            vcf_status = variety_info.get("VCF Status")
            print(f"variety_info: {variety_info}")
            fixed_vcf_item_data = {
                'processed_name': processed_name,
                'variety_id': variety_info.get('variety_id') or variety_id,
                'status': vcf_status or 'No VCF'
            }

            # GWAS 탭 on/off
            if vcf_status and vcf_status != 'No VCF':
                initial_vcf_values = [vcf_status]
                initial_analysis_tabs = [
                    dbc.Tab(label="Phenotype Analysis", tab_id="phenotype-tab"),
                    dbc.Tab(label=f"GWAS Analysis", tab_id="gwas-tab", disabled=False),
                ]
            else:
                initial_analysis_tabs = [
                    dbc.Tab(label="Phenotype Analysis", tab_id="phenotype-tab"),
                ]

            # ✅ 통합 패널로 phenotype 초기 콘텐츠 생성 (PHENO_DF 기반)
            safe_it = str(
                variety_info.get("IT Number")
                or variety_info.get("it_number")
                or variety_info.get("variety_id")
                or processed_name
            )
            print(f"safe_it: {safe_it}")
            if PHENO_DF is None or PHENO_DF.empty:
                initial_tab_content = html.Div([
                    dbc.Alert([
                        html.H5("❌ No Phenotype Data", className="alert-heading"),
                        html.P("PHENO_DF is empty or not loaded.")
                    ], color="danger")
                ])
            else:
                #print(f"safe_it: {safe_it}")
                initial_tab_content = create_phenotype_panel_unified(
                    selected_nodes=[safe_it],
                    pheno_df=PHENO_DF,
                    table_enabled=True
                )

            initial_active_tab = "phenotype-tab"

        except Exception as e:
            print(f"❌ Auto-initialization error: {e}")
            initial_tab_content = html.Div([
                dbc.Alert([
                    html.H5("Phenotype Analysis Ready", className="alert-heading"),
                    html.P("Search for a variety and run analysis.")
                ], color="success")
            ])
    # 4. 완전한 레이아웃 생성 (좌우 분할 레이아웃)
    return html.Div([
            
             
        
        # CSS for Discovery span hover effects - improved formatting
        html.Link(
            rel='stylesheet',
            href='data:text/css,' + 
                 '.ecotype-span-hover:hover { ' +
                 'background-color: #007bff !important; ' +
                 'color: white !important; ' +
                 'border-color: #007bff !important; ' +
                 'transform: scale(1.05); ' +
                 'box-shadow: 0 2px 8px rgba(0,123,255,0.3); ' +
                 '} ' +
                 '.variety-group-span-hover:hover { ' +
                 'background-color: #28a745 !important; ' +
                 'color: white !important; ' +
                 'border-color: #28a745 !important; ' +
                 'transform: scale(1.05); ' +
                 'box-shadow: 0 2px 8px rgba(40,167,69,0.3); ' +
                 '}'
        ),
        
        # Header
        create_header(),
        
        # 📍 품종 정보 카드 영역 (초기값 설정)
        html.Div(id='header-section', children=header_section),
        
        # URL 알림 영역 (초기값 설정)
        html.Div(id='url-notification', children=html.Div(), style={'margin': '10px'}),

       
        html.Div([

            
            # 좌측 패널: 계보도 시각화 (40% 너비)
            html.Div([
                # Pedigree 헤더와 접기 버튼
                html.Div(
                    [
                        html.H5(
                            "🌳 Pedigree Visualization",
                            id="pedigree-title",
                            style={
                                "margin": "0",
                                "fontWeight": "bold",
                                "flex": "1 1 auto",
                                "whiteSpace": "nowrap",
                                "overflow": "hidden",
                                "textOverflow": "ellipsis"
                            }
                        ),
                        dbc.Button(
                            html.I(className="fas fa-chevron-left"),
                            id="pedigree-collapse-button",
                            size="sm",
                            color="primary",
                            outline=True,
                            style={
                                "flex": "0 0 auto",
                                "border": "2px solid #007bff",
                                "borderRadius": "6px",
                                "padding": "4px 8px",
                                "transition": "all 0.2s ease"
                            },
                            title="collapse/expand"
                        )
                    ],
                    id="pedigree-header",
                    style={
                        "display": "flex",
                        "alignItems": "center",
                        "justifyContent": "space-between",
                        "gap": "8px",
                        "marginBottom": "15px",
                        "flexWrap": "nowrap"
                    }
                ),
                
                # 접기 가능한 내용 영역
                html.Div([
                
                # 제어 패널
                html.Div([
                    # 모드 표시 (고정)
                    

                    # 검색 (nopedi case에서만 표시) - 자동완성 기능 포함 - 양옆 배치
                    html.Div([

                            html.Div([
                                html.P(
                                    "Discover pedigree connections by adding varieties with pedigree data.",
                                    style={'fontWeight': '500', 'fontSize': '0.9rem', 'marginBottom': '4px'}
                                ),
                                html.P([
                                    "Use the ",
                                    html.B("‘Add Additional Varieties’"),
                                    " button below to explore varieties with available pedigree data."
                                ],
                                style={
                                    'color': '#555',
                                    'fontSize': '0.85rem',
                                    'marginBottom': '6px'
                                }),
                                html.Div(id="addv-pedi-message", className="mt-2")
                            ],
                            style={
                                'backgroundColor': '#f8f9fa',
                                'border': '1px solid #dee2e6',
                                'borderRadius': '6px',
                                'padding': '10px',
                                'marginTop': '10px'
                            })
                    ], 
                    id="nopedi-guide-container",
                    className="mb-3", style={'display': 'block' if processed_name == 'nopedi' else 'none'}),

                    # 🚩 Reset View 버튼을 pedigree 캔버스로 이동하여 제거


                    # Unified Pedigree Path (works for both nopedi and regular cases)
                    html.Div([
                        html.H6("Pedigree Set", className="mb-2"),
                        html.Div(id='pedigree-breadcrumb', 
                                children=[],
                                style={'minHeight': '30px'})
                    ], id='pedigree-path-container', style={
                        'padding': '8px 10px',
                        'backgroundColor': '#e8f4fd',
                        'borderRadius': '6px',
                        'border': '1px solid #b8daff',
                        'borderLeft': '3px solid #007bff',
                        'marginBottom': '10px',
                        'display': 'block'  # Always visible for unified approach
                    }),
                html.Div(
                    [html.Div(additional_varieties_modal)],
                    style={
                        'display': 'flex',
                        'justifyContent': 'flex-end',  # ✅ 오른쪽 정렬
                        'marginTop': '8px'
                    }
                ),
                    # Unified Selected Varieties (Apply 전) - works for both nopedi and regular cases
                   
                    
                    

                ], style={
                    'padding': '15px',
                    'backgroundColor': '#ffffff',
                    'borderRadius': '8px',
                    'border': '1px solid #dee2e6',
                    'marginBottom': '15px'
                }),
                
                html.Div(
                    [
                        # 🎨 Phenotype Legend
                        html.Div(
                            [
                                html.Div("🎨 Phenotype Legend", style={
                                    "fontWeight": "bold",
                                    "fontSize": "13px",
                                    "marginBottom": "4px",
                                    "color": "#343a40",
                                }),
                                html.Div(
                                    id="trait-legend-container",
                                    style={
                                        "display": "flex",
                                        "flexDirection": "column",
                                        "gap": "4px",
                                        "maxHeight": "160px",
                                        "overflowY": "auto",
                                        "paddingRight": "6px",
                                        "width": "100%",
                                    }
                                )
                            ],
                            id="trait-legend-section",
                            style={
                                "display": "none",   # 기본은 숨김
                                "padding": "15px",
                                "backgroundColor": "#ffffff",
                                "borderRadius": "8px",
                                "border": "1px solid #dee2e6",
                                "marginBottom": "15px",
                            },
                        ),

                        # 🧬 Genotype Legend (GWAS)
                        html.Div(
                            [
                                html.Div("🧬 Genotype Legend", style={
                                    "fontWeight": "bold",
                                    "fontSize": "13px",
                                    "marginBottom": "4px",
                                    "color": "#343a40",
                                }),
                                html.Div(
                                    id="gwas-genotype-legend",
                                    style={
                                        "display": "flex",
                                        "flexDirection": "column",
                                        "gap": "4px",
                                        "maxHeight": "160px",
                                        "overflowY": "auto",
                                        "paddingRight": "6px",
                                        "width": "100%",
                                    }
                                )
                            ],
                            id="gwas-genotype-legend-section",
                            style={
                                "display": "none",   # 기본 숨김
                                "padding": "15px",
                                "backgroundColor": "#ffffff",
                                "borderRadius": "8px",
                                "border": "1px solid #dee2e6",
                                "marginBottom": "15px",
                            },
                        ),
                    ],
                    style={
                        "width": "100%",  # ✅ 같은 부모 내부에 배치
                    }
                ),
                # 계보도 시각화 (상대적 위치 컨테이너)
                html.Div([
                    # 🚩 Reset Pedigree 버튼 (깔끔한 아이콘 버튼)
                    html.Div(
                        [
                            # 상단 바 (Legend + Buttons)
                            html.Div(
                                [
                                    # 🔹 좌측 Legend
                                    html.Div(
                                        [
                                            html.H6(
                                                "🏷️ Pedigree Legend",
                                                style={
                                                    "margin": "0 10px 0 0",
                                                    "fontSize": "12px",
                                                    "fontWeight": "bold",
                                                    "color": "#f8f9fa",
                                                },
                                            ),
                                            _legend_item("#1C6BA0", "circle", "Genotype"),
                                            _legend_item("#27ae60", "circle", "Phenotype"),
                                            _legend_item("#e67e22", "circle", "Phenotype + Genotype"),
                                            _legend_item("#bdc3c7", "circle", "No Data"),
                                            _legend_ring_item("#ff33cc", "circle", "Selected"),  # 🌿 외곽만 초록색 #2ecc71->ff33cc
                                            _legend_item("transparent", "circle", "In Pedi", border="#f8f9fa"),
                                            _legend_item("transparent", "rectangle", "In Non-Pedi", border="#f8f9fa"),
                                        ],
                                        style={
                                            "display": "flex",
                                            "alignItems": "center",
                                            "gap": "4px",
                                        },
                                    ),

                                    # 🔹 우측 버튼 4개
                                    html.Div(
                                        [
                                            html.Button(html.I(className="fas fa-expand-arrows-alt"),
                                                        id="btn-expand", title="Expand node(s)",
                                                        n_clicks=0,
                                                        style={"background": "#343a40", "color": "white",
                                                            "border": "none", "padding": "6px 8px",
                                                            "borderRadius": "4px", "cursor": "pointer"}),
                                            html.Button(html.I(className="fas fa-trash"),
                                                        id="btn-remove", title="Remove node(s)",
                                                        n_clicks=0,
                                                        style={"background": "#343a40", "color": "white",
                                                            "border": "none", "padding": "6px 8px",
                                                            "borderRadius": "4px", "cursor": "pointer"}),
                                            html.Button(html.I(className="fas fa-sync-alt"),
                                                        id="btn-reset", title="Reset View",
                                                        n_clicks=0,
                                                        style={"background": "#343a40", "color": "white",
                                                            "border": "none", "padding": "6px 8px",
                                                            "borderRadius": "4px", "cursor": "pointer"}),
                                            html.Button(html.I(className="fas fa-expand"),
                                                        id="btn-wide", title="Wide view",
                                                        n_clicks=0,
                                                        style={"background": "#343a40", "color": "white",
                                                            "border": "none", "padding": "6px 8px",
                                                            "borderRadius": "4px", "cursor": "pointer"}),
                                        ],
                                        style={
                                            "display": "flex",
                                            "alignItems": "center",
                                            "gap": "6px",
                                        },
                                    ),
                                ],
                                style={
                                    "display": "flex",
                                    "justifyContent": "space-between",
                                    "alignItems": "center",
                                    "backgroundColor": "rgba(52,58,64,0.9)",
                                    "padding": "6px 12px",
                                    "borderRadius": "8px",
                                    "boxShadow": "0 2px 6px rgba(0,0,0,0.25)",
                                },
                            ),

                            html.Div(
                                [
                                    # 🔹 상단: component-badge-container + Data Info 버튼 (한 줄 정렬)
                                    html.Div(
                                        [ html.Div(),
                                            # 왼쪽: Pedigree component badges
                                            html.Div(
                                                id="component-badge-container",
                                                style={
                                                    "display": "flex",
                                                    "justifyContent": "center",  # ✅ 중앙 정렬
                                                    "flexWrap": "wrap",
                                                    "gap": "6px",
                                                    "alignItems": "center",
                                                    "flex": "2",  # 중간 영역을 넓게
                                                },
                                            ),

                                            # 오른쪽: Data Information 버튼
                                            html.Button(
                                                [
                                                    html.I(className="fas fa-database", style={"marginRight": "6px"}),  # 💾 아이콘
                                                    "Data Information"
                                                ],
                                                id="btn-profiles-toggle",
                                                title="Show detailed profiles",
                                                n_clicks=0,
                                                style={
                                                    "display": "flex",
                                                    "alignItems": "center",
                                                    "gap": "6px",
                                                    "backgroundColor": "#007bff",
                                                    "color": "white",
                                                    "border": "none",
                                                    "borderRadius": "6px",
                                                    "padding": "4px 10px",
                                                    "cursor": "pointer",
                                                    "fontSize": "13px",
                                                    "fontWeight": "500",
                                                    "boxShadow": "0 1px 2px rgba(0,0,0,0.15)",
                                                    "transition": "all 0.2s ease-in-out",
                                                    #"flex": "1",  # 버튼은 오른쪽 영역
                                                    "justifyContent": "flex-end",
                                                },
                                            ),
                                        ],
                                        style={
                                             "display": "flex",
                                            "alignItems": "center",
                                            "justifyContent": "space-between",
                                            "width": "100%",
                                            "marginTop": "6px",
                                        },
                                    ),
                                    # 🔹 Trait Color Mapping Panel (일반 패널 형태)
                                    

                                    # 🔹 아래: profiles panel (토글)
                                    html.Div(
                                        id="profiles-panel",
                                        style={
                                            "display": "none",
                                            "backgroundColor": "#f8f9fa",
                                            "border": "1px solid #dee2e6",
                                            "borderRadius": "6px",
                                            "marginTop": "6px",
                                            "padding": "10px",
                                            "maxHeight": "220px",
                                            "overflowY": "auto",
                                            "boxShadow": "inset 0 1px 3px rgba(0,0,0,0.1)",
                                        },
                                        children=[
                                            html.H6("📋 Detailed Profiles", style={
                                                "fontWeight": "bold",
                                                "color": "#495057",
                                                "marginBottom": "8px"
                                            }),
                                            html.Div(id="profiles-content", children=[
                                                html.P("No profiles loaded.", style={"color": "#6c757d", "fontStyle": "italic"})
                                            ]),
                                        ],
                                    ),
                                ]
                            ),


                        ],
                        style={
                            "position": "absolute",
                            "top": "10px",
                            "left": "10px",
                            "right": "10px",
                            "zIndex": 1000,
                        },
                    ),
                        
                            html.Div(id="wide-overlay", style={
                                'position': 'fixed',
                                'top': 0, 'left': 0,
                                'width': '100%', 'height': '100%',
                                'backgroundColor': 'rgba(0,0,0,0.5)',
                                'zIndex': 999,
                                'display': 'none'
                            }),
                    
                    
                   

                    cyto.Cytoscape(
                        id='pedigree-cytoscape',
                        style={'height': '600px', 'width': '100%'},
                        elements=initial_pedigree_elements,
                        stylesheet=get_default_stylesheet2(processed_name),  # ✅ 기본 스타일 복구
                        layout={
                            'name': 'dagre',
                            'rankDir': 'TB',
                            'ranker': 'network-simplex',
                            'animate': True,
                            'animationDuration': 500,
                            'fit': True,
                            'padding': 75,
                            'spacingFactor': 1.2,
                            'nodeDimensionsIncludeLabels': True,
                            'rankSep': 120,
                            'nodeSep': 100,
                            'edgeSep': 50,
                            'avoidOverlap': True

                        },
                        userZoomingEnabled=True,
                        userPanningEnabled=True,
                        boxSelectionEnabled=True,
                        autoungrabify=False,
                        autounselectify=False
                    ), 
                    html.Div(
                        id='tooltip',
                        style={
        'position': 'absolute',
        'bottom': '300px',
        'left': '10px',
        'zIndex': 1000,
        'backgroundColor': 'rgba(0, 0, 0, 0.85)',
        'color': 'white',
        'border': '1px solid #333',
        'borderRadius': '6px',
        'padding': '8px 10px',
        'boxShadow': '0 2px 8px rgba(0,0,0,0.3)',
        'fontSize': '12px',
        'maxWidth': '260px',
        'display': 'none',
        'visibility': 'hidden',
        'opacity': '0',
        'transition': 'opacity 0.15s ease-in-out',
        'pointerEvents': 'none',  # 툴팁이 마우스 이벤트 가로채지 않도록
        'willChange': 'opacity, transform',
    }
                    ),
                       dcc.Store(id='tooltip-sync'),


#tooltip 위치
                    



                    # 🚩 간소화된 JavaScript - tooltip은 Python callback으로 처리됨

                ], style={
                    'padding': '15px',
                    'backgroundColor': '#ffffff',
                    'borderRadius': '8px',
                    'border': '1px solid #dee2e6',
                    'boxShadow': '0 2px 4px rgba(0,0,0,0.05)',
                    'position': 'relative',  # sticky에서 relative로 변경하여 tooltip positioning 지원
                   # 'top': '20px',  # 상단에서 20px 떨어진 위치에 고정
                    'backgroundColor': 'white', # 배경 덮어쓰기
                    'zIndex': 10      

                }, id='pedigree-main-container')
                
                ], id="pedigree-collapse", style={"display": "block",'position':'sticky','top':'25px'})  # 접기 가능한 내용 영역 닫기

,
               


            ], id={'type': 'left-panel', 'index': 0}, style={
                'width': '40%',
                'paddingRight': '10px',
                'overflow': 'visible',  # ✅ 반드시 visible
                'position': 'relative', # 기본값 유지
            }),
            
            # 우측 패널: 분석 (60% 너비)
            html.Div([
                html.Div(
                    [
                        html.H5(
                            "📊 Analysis Dashboard",
                            id="analysis-dashboard-title",
                            style={
                                "margin": "0",
                                "fontWeight": "bold",
                                "flex": "1 1 auto",
                                "whiteSpace": "nowrap",
                                "overflow": "hidden",
                                "textOverflow": "ellipsis",
                            },
                        )
                    ],
                    id="analysis-dashboard-header",
                    style={
                        "display": "flex",
                        "alignItems": "center",
                        "justifyContent": "space-between",
                        "gap": "8px",
                        "marginBottom": "15px",
                        "flexWrap": "nowrap",
                    },
                ),
                
                html.Div([
                    # Tab system
                    dbc.Tabs(
                        id="analysis-tabs",
                        active_tab=initial_active_tab,
                        children=initial_analysis_tabs
                    ),
                    
                    html.Div(id="tab-content", children=initial_tab_content, style={'margin-top': '20px'})
                ], style={
                    'padding': '15px',
                    'backgroundColor': '#ffffff',
                    'borderRadius': '8px',
                    'border': '1px solid #dee2e6',
                    'boxShadow': '0 2px 4px rgba(0,0,0,0.05)'
                })
                
            ], id={'type': 'right-panel', 'index': 0}, style={
                'width': '60%',
                'paddingLeft': '10px',
                'overflow': 'auto'
            })
            
        ], style={
            'display': 'flex',
            'margin': '20px',
            'gap': '20px'
        }),
        
        dcc.Store(id='gwas-selected-label-store', data={}),
        
        dcc.Store(id='phenotype-label-store', data={}),
        dcc.Store(id='component-info-store', data=[]),
        dcc.Store(id='highlighted-components-store', data=[]),
        dcc.Store(id='selected-nodes-store', data=[]),
        
        

        dcc.Store(id='phenotype-nopedi-data-store', data=None),
        dcc.Store(id='pedigree-path-child-store',data={}),
        # Store components
        dcc.Store(id='snp-occurrence-store', data={"enabled": False, "threshold": None}),
        dcc.Store(id='filter-meta-store', data={}),
        dcc.Store(id='filter-meta-store2', data={}),
        #dbc.Collapse(id='selected-traits-collapse', is_open=False),
        #dcc.Store(id="selected-traits-collapse-store", data=False),
        
        dcc.Store(id='maf-filter-store', data={'enabled': False, 'cutoff': 0.05}),
        dcc.Store(id='selected-node-store', data=None),
        dcc.Store(id='multi-selected-nodes', data=[]),  # 호환성을 위해 유지 (빈 배열)
        dcc.Store(id='available-nodes-store', data=[]),  # elements 스캔 결과(가용 노드) 저장
        dcc.Store(id='selection-mode', data=False),
        dcc.Store(id='current-vcf-values', data=initial_vcf_values),
        dcc.Store(id='gwas-selected-traits-store', data=[]),
        dcc.Store(id='gwas-filter-states', data={}),
        dcc.Store(id='gwas-input-loci-store', data=[]),
        
        dcc.Store(id='pedigree-path-store', data=[]),  # pedigree 확장 경로 저장소
        dcc.Store(id='trait-color-store', data={}),     # 트레이트 컬러 저장소 추가
        dcc.Store(id='fixed-vcf-item', data=fixed_vcf_item_data or {}),
        dcc.Store(id='fixed-vcf-item-nopedi', data=fixed_vcf_item_nopedi_data or {}),
        dcc.Store(id='search-selected-variety-name', data=None),  # 사용자가 검색/확정한 품종명
        dcc.Store(id='pedigree-elements-store', data=[]),  # 검색한 이름 기준으로 생성한 pedigree elements
        # Integrated Filter System Stores
        dcc.Store(id='integrated-selected-traits-store', data=[]),
        dcc.Store(id='pedigree-click-tracker', data={'last_click': 0, 'node_id': None}),  # 클릭/더블클릭 구분용
        dcc.Store(id='integrated-trait-groups-store', data={'groups': {}, 'max_groups': 9}),
        dcc.Store(id='integrated-secondary-filter-enabled', data=False),
        dcc.Store(id='integrated-active-filter-tab', data=None),
        dcc.Store(id='integrated-filter-b-state', data={'enabled': True, 'pvalue_cutoff': 5}),
        dcc.Store(id='integrated-filter-c-state', data={'enabled': False, 'maf_cutoff': 0.05}),
        dcc.Store(id='integrated-filter-d-state', data={'enabled': False}),
        # GWAS Legend Filter State Store (ref 코드 로직 적용)
        dcc.Store(id='gwas-legend-filters', data={'sample': {}, 'subtrait': {}}),
        
        html.Div(id='subtrait-plot-container',style={'display': 'none'}),
        html.Div(id='subtrait-plot-container-group',style={'display': 'none'}),
        
        # 🚩 Enhanced stores for expand history tracking
       
        dcc.Store(id='child-nodes-store', data=[]),     # 확장으로 추가된 자식 노드들
        
        
        # 🚩 Reset View trigger store
        dcc.Store(id='reset-view-trigger', data=0),
        
        # Note: nopedi auto-detected varieties now unified with available-nodes-store
        
        
        dcc.Store(id="last-hover-ts", data=0.0),
        # 🏗️ GWAS 상태 관리
        dcc.Store(id='gwas-state', data={'status': 'idle', 'used_samples': [], 'last_error': ''}),
        
        # 🖱️ 마우스 버튼 감지를 위한 store
        dcc.Store(id='mouse-button-store', data={'button': 'left', 'timestamp': 0}),
        
        
        # 🎯 Select/Apply 분리를 위한 stores
        dcc.Store(id='applied-varieties-store', data=[]),  # Apply된 품종들 (unified for both regular and nopedi cases)
        #dcc.Store(id='scatter-refresh-trigger', data=0),
        #dcc.Store(id='scatter-refresh-trigger-group', data=0),

        # 🚩 Hidden GWAS components (for callback compatibility in nopedi case)
        #dcc.Store(id='active-samples-store', data=[]),
        dcc.Store(id="gwas-click-store", data=[]),
        dcc.Store(id='gwas-sample-scatter-group-data'),
        dcc.Store(id="pheno-store", data=PHENO_DF.to_dict("records")),
        dcc.Store(id="sample-mode-store",data='each'),
        dcc.Store(id="copy2selected-traits-store",data=[]),
        # Debug div for available-nodes-store (hidden)
        #html.Div(id='debug-available-nodes', style={'display': 'none'}),
        dcc.Store(id="presence-sets-store",data=None),
        dcc.Store(id="available-samples-store"),
        # Hidden stores and tables for callback compatibility
        dcc.Store(id='phenotype-data-store', data=[{'id': variety_id}]),  
        dcc.Store(id='copy2pheno-selected-samples-store', data=[{'id': variety_id}]),
        dcc.Store(
            id='copy2gwas-selected-samples-store',
            data=[{
                'vcf_status': vcf_status,
                'id': parsed_id,  # ✅ 이 부분만 사용!
            }] if vcf_status and vcf_status != '-' else []
        ),
        dcc.Store(id='gwas-selected-samples-store', data=[vcf_status] if vcf_status and vcf_status != '-' else []),
        dcc.Store(id='gwas-combo-store', storage_type='memory'),
        dcc.Store(id="df-vp-store"),
        dcc.Store(id="df-vp-store2"),
        dcc.Store(id='active-trait-store', data={'trait': None, 'category': None}),

        #dcc.Store(id="sampleorder-store"),
        dcc.Store(id="processed-df-store", data={}),   # dict 저장용
        dcc.Store(id="final-df-store", data=None),     # 단일 DataFrame 저장용
        html.Div([       
            dcc.Graph(id='gwas-sample-scatter', figure={}),
        dcc.Graph(id='gwas-sample-scatter-group', figure={}),
             

            dcc.Input(
                                                    id="pvalue-cutoff",
                                                    type="number",
                                                    placeholder="5",
                                                    value=5,
                                                     style={"display": "none"}
                                                ),
            html.Div(id="group-info-display", style={"display": "none"}),
            dcc.Graph(id="integrated-trait-barplot", figure={}),  # 빈 그래프
            #dcc.Graph(id='integrated-trait-barplot-group',  figure={}),
    
            html.Div(id="individual-sample-size", style={"display": "none"}),
            #dbc.Switch(id="group-enable-switch", label="Enable Group Mode", value=False,            style={"display": "none"}),
            dbc.Switch(id='maf-enabled', label="Enable MAF Filter", value=False,
            style={"display": "none"}),
            
            dcc.RadioItems(
            id="sample-mode-radio",
            options=[],
            value=None,
            style={"display": "none"},)
            ,
            html.Div(id="trait-group-container", style={"display": "none"},),
            html.Div(id="variant-group-container", style={"display": "none"},),
            html.Div(id="sample-group-container", style={"display": "none"},),
            html.Div(id="selected-traits-section", style={"display": "none"},),
            # 👇 placeholder들 미리 선언    
            html.Div(id="integrated-selected-traits-table", style={"display": "none"},),
            html.Div(id="selected-traits-toggle", style={"display": "none"},),
            html.Div(id="selected-traits-status" ,style={"display": "none"},),
            html.Pre(id="selected-traits-summary" ,style={"disptoggle_tablelay": "none"},),
       
    #html.Div(id="individual-sample-size", style={"display": "none"}),
    html.Div(id="merge-sample-size", style={"display": "none"}),
    html.Div(id="unique-sample-size", style={"display": "none"}),
    html.Div(id="individual-sample-section", style={"display": "none"}),
        html.Div(id="merge-sample-section", style={"display": "none"}),
        html.Div(id="unique-sample-section", style={"display": "none"}),
        dash_table.DataTable(
            id='mannequin-combined-table-copy2',
            columns=[],
            data=[]
        ),
        dcc.Dropdown(
            id='integrated-subtrait-dropdown',
            options=[],        # 나중에 콜백에서 채움
            value=None,
            placeholder="Select subtraits..."
        ,style={'display': 'none'}),
        dcc.RadioItems(
                                    id='individual-sample-radio-visible',
                                    options=[],  # 동적 업데이트
                                    value=None,
                                    inline=True
                                ),
        dcc.RadioItems(
                                id='group-type-radio',
                                options=[
                                    
                                ],
                                value=None,
                                inline=True,
                                #labelStyle={'margin-right': '10px', 'margin-left': '3px'},
                                #style={'marginBottom': '5px','marginLeft':'15px'}
                            ),

        ], style={"display": "none"})

,
        create_footer()

        # Hidden container for GWAS pedi samples (callback compatibility)
        
    ])

def get_marker_symbol_text(symbol):
    """Plotly 마커 심볼을 텍스트 문자로 변환"""
    symbol_map = {
        'circle': '●',
        'square': '■', 
        'diamond': '♦',
        'cross': '+',
        'x': '×',
        'triangle-up': '▲',
        'triangle-down': '▼',
        'triangle-left': '◄',
        'triangle-right': '►',
        'star': '★',
        'hexagram': '✡',
        'pentagon': '⬟'
    }
    return symbol_map.get(symbol, '●')  # 기본값은 원


def create_custom_legend(fig, selected_samples, show_unique_only, subtrait_color_map, sample_marker_map):
    """
    커스텀 범례 생성 - text 기반으로 title 아래에 배치
    - Row 1: 축 정보 (X-axis, Y-axis)
    - Row 2-3: Subtraits (9개) - 5개씩 2줄
    - Row 4+: Samples (show_unique_only가 True일 때)
    """
    import plotly.graph_objects as go
    
    # 범례 설정 - plot area와 분리된 위치에 배치
    legend_x_start = 0.05  # 시작 X 위치
    # 범례를 margin 영역 안에 배치하되 plot과 겹치지 않도록 조정
    legend_y_start = 0.98  # 시작 Y 위치를 더 위로 이동
    
    cols_per_row = 5
    row_height = 0.02  # 각 행의 높이를 줄여서 더 조밀하게 배치
    col_width = 0.18   # 각 열의 너비
    marker_size = 8
    text_size = 9
    
    # 범례 행 개수 계산 (마진 조정용)
    sample_legend_rows = 0 if not (show_unique_only and selected_samples) else (len(selected_samples) + 4) // 5
    total_rows = 1 + 2 + sample_legend_rows  # 축 정보 1행 + Subtrait 2행 + 샘플 행들
    
    # 0. 축 정보 (1행)
    axis_info_y = legend_y_start - row_height
    
    # X-axis 정보
    fig.add_annotation(
        x=legend_x_start,
        y=axis_info_y,
        text="GWAS-related sample SNP<br><b>X-axis:</b> Position (bp)",
        showarrow=False,
        font=dict(size=text_size, color='black'),
        xanchor='left',
        yanchor='middle',
        xref='paper',
        yref='paper'
    )
    
    # Y-axis 정보
    fig.add_annotation(
        x=legend_x_start + 0.5,
        y=axis_info_y,
        text="<b>Y-axis:</b> -log₁₀(P)",
        showarrow=False,
        font=dict(size=text_size, color='black'),
        xanchor='left',
        yanchor='middle',
        xref='paper',
        yref='paper'
    )
    
    # 1. Subtrait 범례 (2줄로 배치)
    for i, subtrait in enumerate(ALL_SUBTRAITS):
        row = i // cols_per_row
        col = i % cols_per_row
        
        x_pos = legend_x_start + col * col_width
        y_pos = legend_y_start - (row + 2) * row_height  # 축 정보 행 다음부터 시작 (축 정보 1행 + 여백 고려)
        
        # 마커 추가 (주석으로 배치)
        fig.add_annotation(
            x=x_pos,
            y=y_pos,
            text="●",  # 원 심볼
            showarrow=False,
            font=dict(
                size=marker_size + 6,  # 크기 조절
                color=subtrait_color_map.get(subtrait, '#888888')
            ),
            xanchor='center',
            yanchor='middle',
            xref='paper',
            yref='paper'
        )
        
        # 텍스트 추가 (짧은 이름으로 표시)
        short_name = subtrait.split('(')[0] if '(' in subtrait else subtrait
        fig.add_annotation(
            x=x_pos + 0.02,
            y=y_pos,
            text=short_name,
            showarrow=False,
            font=dict(size=text_size, color='black'),
            xanchor='left',
            yanchor='middle',
            xref='paper',
            yref='paper'
        )
    
    # 2. Sample 범례 (show_unique_only가 True일 때)
    if show_unique_only and selected_samples:
        start_row = 4  # 축 정보 1행 + 여백 1행 + Subtrait 2줄 다음부터 시작
        
        for i, sample in enumerate(selected_samples):
            row = start_row + (i // cols_per_row)
            col = i % cols_per_row
            
            x_pos = legend_x_start + col * col_width
            y_pos = legend_y_start - row * row_height
            
            # 🔥 수정: symbol은 sample ≥2 AND unique filter에서만 사용
            if len(selected_samples) >= 2 and show_unique_only:
                marker_symbol = sample_marker_map.get(sample, 'circle')
            else:
                marker_symbol = 'circle'  # 기본은 항상 circle
            
            # 마커 추가 (샘플용 심볼)
            marker_text = get_marker_symbol_text(marker_symbol)
            fig.add_annotation(
                x=x_pos,
                y=y_pos,
                text=marker_text,
                showarrow=False,
                font=dict(
                    size=marker_size + 6,
                    color='gray'
                ),
                xanchor='center',
                yanchor='middle',
                xref='paper',
                yref='paper'
            )
            
            # 텍스트 추가
            fig.add_annotation(
                x=x_pos + 0.02,
                y=y_pos,
                text=f"Sample: {sample}",
                showarrow=False,
                font=dict(size=text_size, color='black'),
                xanchor='left',
                yanchor='middle',
                xref='paper',
                yref='paper'
            )




def setup_gwas_callbacks():
    """GWAS 모듈의 모든 콜백들 설정"""
    
    # Pedigree collapse 버튼 콜백 - 사이드바 형태로 수정
    @callback(
        [
            Output('pedigree-collapse', 'style'),
            Output('pedigree-collapse-button', 'children'),
            Output({'type': 'left-panel', 'index': 0}, 'style'),
            Output({'type': 'right-panel', 'index': 0}, 'style'),
            Output('pedigree-title', 'children'),
            Output('pedigree-title', 'style'),
            Output('pedigree-header', 'style'),
        ],
        [Input('pedigree-collapse-button', 'n_clicks')],
        [State('pedigree-collapse', 'style')],
        prevent_initial_call=False
    )
    def toggle_pedigree_collapse(n_clicks, current_style):
        base_style = {
            "position": "sticky",
            "top": "25px",
            "zIndex": 10,
            "background": "white"
        }

        # ✅ 처음 로드 시에는 콜백 실행 금지
        if not n_clicks:
            return (
                {**base_style, "display": "block"},
                html.I(className="fas fa-chevron-left"),
                {'width': '40%', 'paddingRight': '10px'},
                {'width': '60%', 'paddingLeft': '10px'},
                "🌳 Pedigree Visualization",
                {"fontSize": "1.25rem", "fontWeight": "bold", "margin": "0"},
                {
                    "display": "flex",
                    "alignItems": "center",
                    "justifyContent": "space-between",
                    "gap": "8px",
                    "marginBottom": "15px",
                    "flexDirection": "row"
                }
            )

        # 이후부터는 진짜 토글 동작
        is_hidden = current_style and current_style.get("display") == "none"

        if is_hidden:
            # 펼치기
            return (
                {**base_style, "display": "block"},
                html.I(className="fas fa-chevron-left"),
                {'width': '40%', 'paddingRight': '10px'},
                {'width': '60%', 'paddingLeft': '10px'},
                "🌳 Pedigree Visualization",
                {"fontSize": "1.25rem", "fontWeight": "bold", "margin": "0"},
                {
                    "display": "flex",
                    "alignItems": "center",
                    "justifyContent": "space-between",
                    "gap": "8px",
                    "marginBottom": "15px",
                    "flexDirection": "row"
                }
            )
        else:
            # 접기
            return (
                {**base_style, "display": "none"},
                html.I(className="fas fa-chevron-right"),
                {
                    'width': '50px',
                    'paddingRight': '5px',
                    'min-width': '50px',
                    'flex-shrink': '0'
                },
                {
                    'width': 'calc(100% - 60px)',
                    'paddingLeft': '10px',
                    'flex': '1'
                },
                "🌳 pedigree",
                {"fontSize": "0.8rem", "fontWeight": "normal", "margin": "0"},
                {
                    "display": "flex",
                    "alignItems": "center",
                    "justifyContent": "center",
                    "gap": "4px",
                    "marginBottom": "8px",
                    "flexDirection": "column"
                }
            )

        
    @callback(
        [
            Output('pedigree-main-container', 'style'),
            Output('pedigree-path-container', 'style'),
        ],
        [
            Input('fixed-vcf-item-nopedi', 'data'),
            Input('addv-apply-groups-pedi', 'data'),
        ],
        [
            State('pedigree-main-container', 'style'),
            State('pedigree-path-container', 'style'),
        ],
        prevent_initial_call=False
    )
    def control_nopedi_pedigree_main_container(
        nopedi_data, group_pedi,
        main_container_style, pedigree_path_container_style
    ):
        """
        nopedi 상태에서는 group_pedi가 없으면 초기 렌더부터 컨테이너 숨김.
        pedi/normal 상태에서는 항상 표시.
        """
        main_style = dict(main_container_style or {})
        path_style = dict(pedigree_path_container_style or {})

        # ✅ 상태 판별
        is_nopedi_case = bool(
            nopedi_data and isinstance(nopedi_data, dict) and nopedi_data.get("variety_id")
        )
        has_group_pedi = bool(group_pedi and len(group_pedi) > 0)

        print("🧪 nopedi_data:", nopedi_data)
        print("🧪 is_nopedi_case:", is_nopedi_case)
        print("📦 group_pedi:", group_pedi)
        print("✅ has_group_pedi:", has_group_pedi)

        # --- 1️⃣ nopedi case
        if is_nopedi_case:
            if has_group_pedi:
                print("🟩 nopedi + pedi group detected → 컨테이너 표시")
                return (
                    {**main_style, 'display': 'block'},
                    {**path_style, 'display': 'block'},
                )
            else:
                print("🟥 nopedi 초기 상태 → 컨테이너 숨김 (group 없음)")
                return (
                    {**main_style, 'display': 'none'},
                    {**path_style, 'display': 'none'},
                )

        # --- 2️⃣ 일반 (pedi or normal) case
        print("🟦 일반 pedi case → 항상 표시")
        return (
            {**main_style, 'display': 'block'},
            {**path_style, 'display': 'block'},
        )


    # Integrated Filter 메인 토글 콜백
    
    
    # Sample Summary 컨텐츠 업데이트 콜백들
    '''
    @callback(
        Output('sample-summary-variant-profile-content', 'children',allow_duplicate=True),
        [Input('sample-summary-variant-profile-collapse', 'is_open'),
         Input('mannequin-combined-table', 'data'),
        Input('fixed-vcf-item', 'data'),
         Input('fixed-vcf-item-nopedi', 'data')],
        prevent_initial_call=True
    )
    def update_sample_summary_variant_profile_content(is_open, combined_data, fixed_vcf_data, fixed_vcf_nopedi_data):
        """Sample Summary (variant profile) 컨텐츠 업데이트"""
        if not is_open:
            return []
        
        # 샘플 목록 추출
        all_samples = []
        
        # combined_data에서 샘플 목록 추출
        if combined_data:
            for row in combined_data:
                entry_no = row.get('Entry_No.')
                if entry_no:
                    all_samples.append(str(entry_no))
        
        # 초기 데이터가 없으면 fixed_vcf_value 사용
        if not all_samples:
            if fixed_vcf_data and isinstance(fixed_vcf_data, dict):
                fixed_vcf_value = fixed_vcf_data.get('status')
                if fixed_vcf_value:
                    all_samples.append(str(fixed_vcf_value))
            elif fixed_vcf_nopedi_data and isinstance(fixed_vcf_nopedi_data, dict):
                fixed_vcf_value = fixed_vcf_nopedi_data.get('status')
                if fixed_vcf_value:
                    all_samples.append(str(fixed_vcf_value))
        
        all_samples = list(set(all_samples))
        
        if not all_samples:
            return [html.P("No samples available", className="text-muted text-center p-3")]
        
        # Summary 데이터 생성
        summary_df = create_sample_summary_data(all_samples)
        if summary_df is None or summary_df.empty:
            return [html.P("No summary data available for selected samples", 
                          className="text-muted text-center p-3")]
        
        return [
                html.Div([
                    # 추가샘플과함께분석 span
                
                    dbc.Button(
                        [
                            html.I(className="fas fa-plus-circle", style={'marginRight': '6px'}),
                            "Analyze with Additional Samples"
                        ],
                        id='multi-sample-toggle-btn',
                        n_clicks=0,
                        color="secondary",
                        outline=True,
                        style={
                            'fontWeight': '600',
                            'fontSize': '0.9rem',
                            'borderRadius': '8px',
                            'padding': '6px 12px',
                            'marginBottom': '5px'
                        }
                    ),
                    html.Small(
                        f"(Currently {len(all_samples)} samples)",
                        style={'color': '#6c757d', 'fontSize': '0.8rem', 'marginLeft': '6px'}
                    ),
                
                # 테이블/바플롯 전환 버튼
                dbc.ButtonGroup([
                    dbc.Button("Table View", id="variant-profile-table-btn", color="primary", size="sm"),
                    dbc.Button("Bar Plot", id="variant-profile-barplot-btn", color="outline-primary", size="sm")
                ], style={
                    'margin-bottom': '15px',
                    'display': 'flex',
                    'justifyContent': 'center'
                }),
                
                # 컨텐츠 영역
                html.Div([
                    create_sample_summary_table(summary_df, all_samples)
                ], id='variant-profile-content', style={
                    'backgroundColor': '#ffffff',
                    'border': '1px solid #e9ecef',
                    'borderRadius': '4px',
                    'padding': '10px',
                    'boxShadow': '0 1px 3px rgba(0,0,0,0.1)'
                })
                ,html.Div(id='trait-variant-detail-container', style={'marginTop': '15px'})
          
          
          
            ], style={
                'width': '100%',
                'maxWidth': '95%',
                'margin': '0 auto'
            })
        ]

    @callback(
        Output('sample-summary-snp-occurrence-content', 'children'),
        [
            Input('sample-summary-snp-occurrence-collapse', 'is_open'),
            Input('mannequin-combined-table', 'data'),
            Input('fixed-vcf-item', 'data'),
            Input('fixed-vcf-item-nopedi', 'data'),
            Input('snp-occurrence-min-count-store', 'data'),   # ✅ Store를 Input으로
        ],
        prevent_initial_call=False
    )
    def update_sample_summary_snp_occurrence_content(is_open, combined_data, fixed_vcf_data, fixed_vcf_nopedi_data, min_store):
        if not is_open:
            return []

        # --- 샘플 구성 (동일) ---
        all_samples, primary_sample = [], None
        if isinstance(fixed_vcf_data, dict):
            primary_sample = fixed_vcf_data.get('status')
        elif isinstance(fixed_vcf_nopedi_data, dict):
            primary_sample = fixed_vcf_nopedi_data.get('status')
        if primary_sample:
            primary_sample = str(primary_sample)

        if combined_data:
            for row in combined_data:
                entry_no = row.get('Entry_No.')
                if entry_no:
                    all_samples.append(str(entry_no))

        if not all_samples and primary_sample:
            all_samples.append(primary_sample)

        all_samples = list(set(all_samples))
        if primary_sample and primary_sample in all_samples:
            all_samples.remove(primary_sample)
            all_samples.insert(0, primary_sample)

        if not all_samples:
            return [html.P("No samples available", className="text-muted text-center p-3")]

        total_n = len(all_samples)

        # ✅ Store에서 읽은 min_count로 집계 (하드코딩 1 제거)
        try:
            min_count = int(min_store) if min_store is not None else 1
        except Exception:
            min_count = 1
        min_count = max(1, min(min_count, total_n))

        df = create_all_traits_variant_presence_table(all_samples, min_present_count=min_count)
        if df is None or df.empty:
            title_txt = f"Variants ≥ {min_count} / {total_n} samples — 0 rows"
            columns, data = [], []
        else:
            columns, data = _build_table_payload(df, all_samples)
            title_txt = f"Variants ≥ {min_count} / {total_n} samples — {len(data)} rows"

        # ✅ dcc.Input의 value도 Store값으로 렌더 (재그림 시 유지)
        return [
            html.Div([
                html.Div([
                    html.H6("GWAS ID Filter (Count-based)", style={'color': '#2c3e50', 'marginBottom': '12px'}),
                    dbc.Row([
                        dbc.Col([
                            html.Label("Min Present Count", style={'fontWeight': 'bold', 'fontSize': '0.9rem'}),
                            dcc.Input(
                                id='snp-occurrence-min-count',
                                type='number', min=1, max=total_n, step=1, value=min_count,  # ⬅️ 여기!
                                style={'width': '100%', 'marginTop': '5px'}
                            ),
                            html.Small(f"of {total_n} samples",
                                    style={'color': '#6c757d', 'fontSize': '0.75rem'})
                        ], width=4),
                        dbc.Col([
                            dbc.Button("Apply", id='snp-occurrence-apply-btn', color='primary', size='sm',
                                    style={'marginRight': '10px', 'marginTop': '24px'}),
                            dbc.Button("Reset", id='snp-occurrence-reset-btn', color='outline-secondary', size='sm',
                                    style={'marginTop': '24px'}),
                        ], width=8, className="d-flex align-items-end"),
                    ])
                ], style={'border': '2px solid #17a2b8','borderRadius': '8px','padding': '12px',
                        'backgroundColor': '#f0f8ff','marginBottom': '16px'}),

                html.Div([
                    html.H6(title_txt, id='snp-occurrence-table-title',
                            style={'color': '#2c3e50', 'marginBottom': '8px'}),
                    dash_table.DataTable(
                        id='snp-occurrence-table',
                        columns=columns,
                        data=data,
                        sort_action='native',
                        filter_action='native',
                        page_size=5,
                        style_table={'overflowX': 'auto'},
                        style_cell={'fontSize': 13, 'padding': '6px', 'textAlign': 'left'},
                    )
                ], id='snp-occurrence-content',
                style={'backgroundColor': '#ffffff','border': '1px solid #e9ecef','borderRadius': '4px',
                        'padding': '10px','boxShadow': '0 1px 3px rgba(0,0,0,0.1)'}),
            ], style={'width': '100%', 'maxWidth': '95%', 'margin': '0 auto'})]

    
   
    def create_trait_variant_presence_table(trait, vcf_samples):
        gwas_result = get_gwas_data_for_samples_optimized(vcf_samples)
        if not gwas_result:
            return html.Div([html.P("No data for selected samples.", className="text-muted p-2")])

        full_df = gwas_result['full_data']
        if full_df is None or full_df.empty:
            return html.Div([html.P("GWAS data is empty.", className="text-muted p-2")])

        meta_df = full_df[full_df['Trait'] == trait]
        if meta_df.empty:
            return html.Div([html.P("No variants for the selected trait.", className="text-muted p-2")])

        trait_info = get_trait_description(trait)
        trait_id = trait_info.get('trait_id', '')
        mapped_terms = trait_info.get('mapped_terms', '')
        description = trait_info.get('description', '')
        subtrait = meta_df['Subtrait'].iloc[0] if len(meta_df) > 0 else 'Unknown'

        # union key (OSAID 우선)
        def pick_union_key(df):
            for c in ['OSAID','OSA_ID','osaid','OSA Id','OSA id','OSA ID']:
                if c in df.columns:
                    return c
            return 'Variant ID'
        union_key = pick_union_key(full_df)

        sample_data = gwas_result.get('sample_data', {})
        variant_union, per_sample_sets = set(), {}

        for s in vcf_samples:
            sdf = sample_data.get(s)
            if sdf is None or sdf.empty:
                per_sample_sets[s] = set()
                continue
            xs = sdf[(sdf['Trait'] == trait)]
            ids = set(xs[union_key].astype(str).tolist())
            variant_union |= ids
            per_sample_sets[s] = ids

        id_map = {}
        if union_key != 'Variant ID' and 'Variant ID' in meta_df.columns:
            ref = meta_df[[union_key, 'Variant ID']].drop_duplicates()
            id_map = dict(zip(ref[union_key].astype(str), ref['Variant ID'].astype(str)))

        rows = []
        for uk in sorted(variant_union):
            display_vid = id_map.get(str(uk), str(uk))
            present_count = 0
            presence_cols = {}
            for s in vcf_samples:
                has = str(uk) in per_sample_sets.get(s, set())
                presence_cols[s] = 'Present' if has else 'Absent'   # EN only
                if has: present_count += 1

            rows.append({
                'Variant ID': display_vid,
                'Trait': trait,
                'Trait ID': trait_id,
                'Mapped Terms': mapped_terms,
                'Description': description,
                'Subtrait': subtrait,
                **presence_cols,
                'Summary': f'{present_count}/{len(vcf_samples)}'     # EN only
            })

        from dash import dash_table
        columns = (
            [{'name': c, 'id': c} for c in ['Variant ID','Trait','Trait ID','Mapped Terms','Description','Subtrait']]
            + [{'name': s, 'id': s} for s in vcf_samples]
            + [{'name': 'Summary', 'id': 'Summary'}]
        )

        return html.Div([
            html.H6(f"Trait '{trait}' variants (n={len(rows)})"),  # EN title
            dash_table.DataTable(
                id='trait-variant-detail-table',
                columns=columns,
                data=rows,
                sort_action='native',
                filter_action='native',
                page_size=10,
                style_table={'overflowX': 'auto'},
                style_cell={'fontSize': 13, 'padding': '6px'},
            )
        ], className="p-2")


    # 5-1) 바 클릭 → 상세
    @callback(
        Output('trait-variant-detail-container', 'children', allow_duplicate=True),
        Input('sample-summary-barplot', 'clickData'),
        State('mannequin-combined-table', 'data'),
        prevent_initial_call=True
    )
    def update_trait_variant_detail_from_bar(clickData, combined_data):
        if not clickData or 'points' not in clickData or not clickData['points']:
            raise dash.exceptions.PreventUpdate
        trait = clickData['points'][0].get('x')
        
        # combined_data에서 샘플 목록 추출
        all_samples = []
        if combined_data:
            for row in combined_data:
                entry_no = row.get('Entry_No.')
                if entry_no:
                    all_samples.append(str(entry_no))
        all_samples = list(set(all_samples))
        
        return create_trait_variant_presence_table(trait, all_samples)

    def _merge_samples_ordered_unique(fixed_sample, additional_samples):
        fx = fixed_sample or []
        ad = additional_samples or []
        if isinstance(fx, str): fx = [fx]
        if isinstance(ad, str): ad = [ad]
        out = []
        for s in fx + ad:
            if s and s not in out:
                out.append(s)
        return out



    @callback(
        Output('trait-variant-detail-container', 'children'),
        Input('sample-summary-table', 'active_cell'),
        State('sample-summary-table', 'data'),
        State('mannequin-combined-table', 'data'),
        prevent_initial_call=True
    )
    def update_trait_variant_detail_from_table(active_cell, table_data, combined_data):
        if not table_data:
            return html.Div([html.P("No table rows.", className="text-muted p-2")])
        if not active_cell:
            return html.Div([html.P("Click a row to see variant details.", className="text-muted p-2")])

        row_idx = active_cell.get('row')
        if row_idx is None or row_idx < 0 or row_idx >= len(table_data):
            return html.Div([html.P("Invalid selection.", className="text-muted p-2")])

        row = table_data[row_idx]
        trait = row.get('Trait') or row.get('trait') or row.get('Trait ') or row.get('TRAIT')
        if not trait:
            return html.Div([html.P("Trait value not found in the selected row.", className="text-danger p-2")])

        # combined_data에서 샘플 목록 추출
        all_samples = []
        if combined_data:
            for row in combined_data:
                entry_no = row.get('Entry_No.')
                if entry_no:
                    all_samples.append(str(entry_no))
        all_samples = list(set(all_samples))
        if not all_samples:
            return html.Div([html.P("No samples selected.", className="text-muted p-2")])

        return create_trait_variant_presence_table(trait, all_samples)
    '''
    '''
    @callback(
        Output('copy2gwas-selected-samples-store', 'data'),
        [Input('mannequin-combined-table-copy2', 'data'),
         Input('fixed-vcf-item', 'data'),
         Input('fixed-vcf-item-nopedi', 'data')],
        prevent_initial_call=False
    )
    def build_selected_samples_store(combined_data, fixed_vcf_data, fixed_vcf_nopedi_data):
        """결과를 copy2gwas-selected-samples-store.data에 저장"""
        from dash import callback_context
        ctx = callback_context
        

        primary = _extract_fixed_vcf(fixed_vcf_data, fixed_vcf_nopedi_data)

        combined_samples = []
        if combined_data:
            for row in combined_data:
                entry_no = row.get('Entry_No.')
                if entry_no is not None:
                    combined_samples.append(entry_no)

        combined_samples = _dedup_preserve_order(combined_samples)

        if primary:
            rest = [s for s in combined_samples if s != primary]
            final_list = [primary] + rest
        else:
            final_list = combined_samples

        if not final_list and primary:
            final_list = [primary]

        return final_list
    '''

    '''
    # Integrated Filter 개별 필터 토글 콜백들
    @callback(
        [Output('integrated-filter-a-collapse', 'is_open'),
         Output('integrated-filter-b-collapse', 'is_open'), 
         Output('integrated-filter-c-collapse', 'is_open'),
         Output('integrated-filter-d-collapse', 'is_open'),
         Output('integrated-filter-btn-a', 'color'),
         Output('integrated-filter-btn-b', 'color'),
         Output('integrated-filter-btn-c', 'color'), 
         Output('integrated-filter-btn-d', 'color')],
        [Input('integrated-filter-btn-a', 'n_clicks'),
         Input('integrated-filter-btn-b', 'n_clicks'),
         Input('integrated-filter-btn-c', 'n_clicks'),
         Input('integrated-filter-btn-d', 'n_clicks')]
    )
    def toggle_individual_filters(btn_a_clicks, btn_b_clicks, btn_c_clicks, btn_d_clicks):
        ctx = callback_context
        if not ctx.triggered:
            return True, False, False, False, "primary", "secondary", "secondary", "secondary"
        
        trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
        
        # 모든 collapse 상태 초기화
        a_open = b_open = c_open = d_open = False
        a_color = b_color = c_color = d_color = "secondary"
        
        # 클릭된 버튼에 따라 해당 collapse만 열기
        if trigger_id == 'integrated-filter-btn-a':
            a_open, a_color = True, "primary"
        elif trigger_id == 'integrated-filter-btn-b':
            b_open, b_color = True, "primary"  
        elif trigger_id == 'integrated-filter-btn-c':
            c_open, c_color = True, "primary"
        elif trigger_id == 'integrated-filter-btn-d':
            d_open, d_color = True, "primary"
        
        return a_open, b_open, c_open, d_open, a_color, b_color, c_color, d_color
    
    # Pedi related 샘플 테이블 업데이트 콜백 (A-B-C-D 구조 적용)
    def updated_update_pedi_samples_table(pedi_related, vcf_info_df):
        """
        기존 update_pedi_samples_table 함수의 개선 버전
        Mannequin 테이블과 호환되는 형식으로 수정
        """
        if not pedi_related or vcf_info_df.empty:
            return [html.P("No pedi related samples.", style={'color': '#999', 'font-style': 'italic'})]
        
        # VCFinfo.csv에서 pedi_related 샘플들과 Entry_No. 매칭
        vcf_matched_df = vcf_info_df[vcf_info_df['Entry_No.'].isin(pedi_related)]
        
        if vcf_matched_df.empty:
            return [html.P("No matching pedigree samples found in VCF data.", style={'color': '#999', 'font-style': 'italic'})]
        
        # path 상의 파일 존재 여부 확인 (개선된 함수 사용)
        vcf_entry_nos = vcf_matched_df['Entry_No.'].tolist()
        existing_samples = updated_check_sample_files_exist(vcf_entry_nos)
        
        if not existing_samples:
            return [html.P("No sample files found in data path.", style={'color': '#999', 'font-style': 'italic'})]
        
        # path에 존재하는 샘플들만 필터링
        filtered_df = vcf_matched_df[vcf_matched_df['Entry_No.'].isin(existing_samples)]
        
        # 체크박스 기능이 포함된 테이블 생성
        return create_mannequin_table_with_checkboxes(filtered_df, 'pedi_samples')
'''
    
    
    
    # 체크박스 상태를 통합 selected samples 리스트로 변환하는 콜백
    
    
    
    
    # B: Pedi related store 업데이트 (expand로 추가된 것들, nopedi case 포함)
    
    
    # C: Discovery available store 업데이트
    
    
    # D: Checked samples store 업데이트 (최종 선택)
    
    
    # 최종 scatter plot용 samples 업데이트 (A + D)
    

    # GWAS 탭 활성화 시 Fixed sample을 selected samples store에 자동 추가
   

    # === 🔥 범례 클릭 상호작용 콜백 (ref 코드 로직 적용) ===
    '''
    @callback(
        Output('gwas-sample-scatter', 'figure', allow_duplicate=True),
        Output('gwas-legend-filters', 'data'),
        Input('gwas-sample-scatter', 'restyleData'),
        State('gwas-legend-filters', 'data'),
        State('copy2gwas-selected-samples-store', 'data'),
        State('integrated-subtrait-dropdown', 'value'),
        State('gwas-selected-traits-store', 'data'),
        State('gwas-filter-states', 'data'),
        State('gwas-input-loci-store', 'data'),
        State('integrated-filter-b-state', 'data'),
        #State('gwas-unique-position-toggle', 'value'),
        prevent_initial_call=True
    )
    def handle_legend_click(restyleData, legend_filters, selected_samples, selected_subtraits, 
                           selected_traits, filter_states, input_loci, pvalue_cutoff):
        """범례 클릭 이벤트를 처리하여 sample/subtrait 가시성 상태 업데이트"""
        # 기본값 보정
        show_unique_only=False
        legend_filters = legend_filters or {'sample': {}, 'subtrait': {}}
        
        # 초기 상태 설정
        if not legend_filters.get('sample') and selected_samples:
            legend_filters['sample'] = {sample: True for sample in selected_samples}
        if not legend_filters.get('subtrait') and selected_subtraits:
            legend_filters['subtrait'] = {subtrait: True for subtrait in (selected_subtraits or [])}
        
        # restyleData가 없으면 현재 상태로 그래프 생성
        if not restyleData or not selected_samples:
            if selected_samples:
                
                gwas_result = get_gwas_data_for_samples_optimized(selected_samples)
                if gwas_result:
                    fig = create_enhanced_gwas_scatter_plot(
                        gwas_result, selected_samples, selected_subtraits, selected_traits,
                        filter_states, input_loci, pvalue_cutoff, show_unique_only,
                        sample_visible=legend_filters.get('sample', {}),
                        subtrait_visible=legend_filters.get('subtrait', {})
                    )
                    return fig, legend_filters
            return no_update, legend_filters
        
        # restyleData 파싱: ([{'visible': ['legendonly']}], [trace_index])
        try:
            changes, indices = restyleData
            new_visible = changes.get('visible', [None])[0]  # True 또는 'legendonly'
            toggled_idx = indices[0]
        except Exception:
            return no_update, legend_filters
        
        # === 🔥 새로운 조건별 범례에 맞는 trace 인덱스 계산 ===
        # 구조: [조합 traces...] + [조건별 더미 traces...]
        
        # 조합 trace 개수 계산
        n_combo = 0
        if selected_samples and selected_subtraits:
            n_combo = len(selected_samples) * len(selected_subtraits) * 12  # chromosome × sample × subtrait
        elif selected_samples:
            n_combo = len(selected_samples) * 12
        else:
            n_combo = 12
            
        current_idx = n_combo
        
        # 조건별 더미 trace 인덱스 계산
        subtrait_dummy_start = current_idx
        n_subtrait = len(selected_subtraits or []) if not (use_group_colors and trait_to_group) else len(trait_to_group or {})
        current_idx += n_subtrait
        
        sample_dummy_start = None
        n_sample = 0
        if len(selected_samples or []) >= 2 and show_unique_only:
            sample_dummy_start = current_idx
            n_sample = len(selected_samples or [])
            current_idx += n_sample
        
        unique_dummy_idx = None
        if show_unique_only:
            unique_dummy_idx = current_idx
            current_idx += 1
            
        maf_dummy_idx = None
        if has_maf_filter:
            maf_dummy_idx = current_idx
            current_idx += 1
        
        print(f"🎯 조건별 Trace 인덱스:")
        print(f"     - 조합: 0-{n_combo-1}")
        print(f"     - subtrait/group: {subtrait_dummy_start}-{subtrait_dummy_start+n_subtrait-1}")
        if sample_dummy_start is not None:
            print(f"     - sample: {sample_dummy_start}-{sample_dummy_start+n_sample-1}")
        if unique_dummy_idx is not None:
            print(f"     - unique: {unique_dummy_idx}")
        if maf_dummy_idx is not None:
            print(f"     - maf: {maf_dummy_idx}")
        
        toggled_on = (new_visible is True)
        
        # 조건별 더미 범례 클릭 처리
        if subtrait_dummy_start <= toggled_idx < subtrait_dummy_start + n_subtrait:
            # Subtrait/Group 더미 범례 클릭
            if use_group_colors and trait_to_group:
                # Group 범례 클릭
                group_names = list(set(group_name for _, (group_name, _) in trait_to_group.items()))
                group_idx = toggled_idx - subtrait_dummy_start
                if 0 <= group_idx < len(group_names):
                    clicked_group = group_names[group_idx]
                    print(f"🎯 Group '{clicked_group}' 클릭: {toggled_on}")
                    # Group 클릭시 해당 group의 모든 trait의 subtrait 조절
                    for trait, (group_name, _) in trait_to_group.items():
                        if group_name == clicked_group:
                            # 해당 trait에 연결된 subtrait들을 찾아서 상태 변경
                            legend_filters['subtrait'][trait] = toggled_on
            else:
                # Subtrait 더미 범례 클릭
                if selected_subtraits:
                    subtrait_idx = toggled_idx - subtrait_dummy_start
                    if 0 <= subtrait_idx < len(selected_subtraits):
                        clicked_subtrait = selected_subtraits[subtrait_idx]
                        legend_filters['subtrait'][clicked_subtrait] = toggled_on
                        print(f"🎯 Subtrait '{clicked_subtrait}' 클릭: {toggled_on}")
                        
                        # 🔥 Unique only filter 사용시 해당 subtrait의 모든 샘플도 함께 제어
                        if show_unique_only and len(selected_samples) >= 2:
                            # 해당 subtrait을 포함하는 모든 샘플의 상태도 함께 변경
                            for sample in selected_samples:
                                if f"{sample}_{clicked_subtrait}" in legend_filters.get('sample_subtrait', {}):
                                    legend_filters['sample_subtrait'][f"{sample}_{clicked_subtrait}"] = toggled_on
                            print(f"🎯 Subtrait '{clicked_subtrait}' 관련 모든 샘플 상태 변경: {toggled_on}")
        
        elif sample_dummy_start is not None and sample_dummy_start <= toggled_idx < sample_dummy_start + n_sample:
            # Sample 더미 범례 클릭 (조건: 2개 이상 + unique filter)
            sample_idx = toggled_idx - sample_dummy_start
            if 0 <= sample_idx < len(selected_samples):
                clicked_sample = selected_samples[sample_idx]
                legend_filters['sample'][clicked_sample] = toggled_on
                print(f"🎯 Sample '{clicked_sample}' 클릭: {toggled_on}")
        
        elif unique_dummy_idx is not None and toggled_idx == unique_dummy_idx:
            # Unique 더미 범례 클릭
            print(f"🎯 Unique 범례 클릭: {toggled_on} (기능 확장 가능)")
            
        elif maf_dummy_idx is not None and toggled_idx == maf_dummy_idx:
            # MAF 더미 범례 클릭
            print(f"🎯 MAF 범례 클릭: {toggled_on} (기능 확장 가능)")
        
        else:
            print(f"🎯 알 수 없는 trace 클릭: index={toggled_idx}")
        
        # 업데이트된 상태로 그래프 다시 생성
        
        gwas_result = get_gwas_data_for_samples_optimized(selected_samples)
        if gwas_result:
            fig = create_enhanced_gwas_scatter_plot(
                gwas_result, selected_samples, selected_subtraits, selected_traits,
                filter_states, input_loci, pvalue_cutoff, show_unique_only,
                sample_visible=legend_filters.get('sample', {}),
                subtrait_visible=legend_filters.get('subtrait', {})
            )
            return fig, legend_filters
        
        return no_update, legend_filters
    '''
    '''
    # Scatter plot 업데이트 (selected samples store 변경에만 반응)
    @callback(
        Output('gwas-sample-scatter', 'figure'),
        Input('copy2gwas-selected-samples-store', 'data'),
        [State('integrated-subtrait-dropdown', 'value'),
         State('gwas-selected-traits-store', 'data'),
         State('gwas-filter-states', 'data'),
         State('gwas-input-loci-store', 'data'),
         State('integrated-filter-b-state', 'data')],
        prevent_initial_call=True
    )
    def update_gwas_scatter(selected_samples, selected_subtraits, selected_traits, filter_states, input_loci, 
                           pvalue_cutoff ):
        """Selected samples 변경 시 scatter plot 업데이트"""
        print(f"🔍 update_gwas_scatter 호출됨: selected_samples={selected_samples} (타입: {type(selected_samples)})")
        show_unique_only=False
        # 초기 로딩 시 빈 figure 반환
        if not selected_samples:
            print(f"❌ 메인 Scatter update skipped: samples={bool(selected_samples)}")
            return {
                'data': [],
                'layout': {
                    'title': '샘플을 선택해주세요',
                    'xaxis': {'title': 'Position'},
                    'yaxis': {'title': '-log10(P-value)'},
                    'annotations': [{
                        'text': '샘플을 선택해 주세요',
                        'xref': 'paper', 'yref': 'paper',
                        'x': 0.5, 'y': 0.5, 'xanchor': 'center', 'yanchor': 'middle',
                        'showarrow': False, 'font': {'size': 16, 'color': 'gray'}
                    }]
                }
            }
        
        print(f"🗺️ 메인 Scatter Plot 업데이트 시작: {selected_samples}")
        print(f"   📊 샘플 개수: {len(selected_samples)}, GWAS 사용 가능: {GWAS_AVAILABLE}")
        
        # selected_samples를 create_enhanced_gwas_scatter_plot이 기대하는 형식으로 변환
        # 함수는 sample ID 문자열 리스트를 기대하므로 변환
        sample_ids = []
        for sample in selected_samples:
            if isinstance(sample, dict):
                sample_ids.append(sample.get('value', sample.get('id', str(sample))))
            else:
                sample_ids.append(str(sample))
        
        print(f"🔍 변환된 샘플 ID: {sample_ids}")
        
        # GWAS 모듈 사용 불가능한 경우
        if not GWAS_AVAILABLE:
            return {
                'data': [],
                'layout': {
                    'title': 'GWAS 모듈 오류',
                    'xaxis': {'title': 'Position'},
                    'yaxis': {'title': '-log10(P-value)'},
                    'annotations': [{
                        'text': 'GWAS 모듈이 사용할 수 없습니다',
                        'xref': 'paper', 'yref': 'paper',
                        'x': 0.5, 'y': 0.5, 'xanchor': 'center', 'yanchor': 'middle',
                        'showarrow': False, 'font': {'size': 16, 'color': 'red'}
                    }]
                }
            }
        
        # GWAS 데이터가 있고 샘플이 있으면 scatter plot 생성 시도
        if sample_ids:
            try:
                print(f"🔒 GWAS scatter plot 생성 시도: {sample_ids}")
                # GWAS 데이터 가져오기
                
                gwas_result = get_gwas_data_for_samples_optimized(sample_ids)
                
                if gwas_result and gwas_result.get('full_data') is not None:
                    return create_enhanced_gwas_scatter_plot(
                        gwas_result, sample_ids, selected_subtraits, selected_traits, 
                        filter_states, input_loci, pvalue_cutoff, show_unique_only,
                        sample_visible={sample: True for sample in sample_ids},
                        subtrait_visible={subtrait: True for subtrait in (selected_subtraits or [])}
                    )
                else:
                    print(f"⚠️ GWAS 데이터를 가져올 수 없음")
                    return {
                        'data': [],
                        'layout': {
                            'title': 'GWAS 데이터 없음',
                            'xaxis': {'title': 'Position'},
                            'yaxis': {'title': '-log10(P-value)'},
                            'annotations': [{
                                'text': '선택된 샘플에 대한 GWAS 데이터가 없습니다',
                                'xref': 'paper', 'yref': 'paper',
                                'x': 0.5, 'y': 0.5, 'xanchor': 'center', 'yanchor': 'middle',
                                'showarrow': False, 'font': {'size': 16, 'color': 'orange'}
                            }]
                        }
                    }
            except Exception as e:
                print(f"❌ Scatter plot 생성 실패: {e}")
                return {
                    'data': [],
                    'layout': {
                        'title': '차트 생성 오류',
                        'xaxis': {'title': 'Position'},
                        'yaxis': {'title': '-log10(P-value)'},
                        'annotations': [{
                            'text': f'Error creating chart: {str(e)}',
                            'xref': 'paper', 'yref': 'paper',
                            'x': 0.5, 'y': 0.5, 'xanchor': 'center', 'yanchor': 'middle',
                            'showarrow': False, 'font': {'size': 14, 'color': 'red'}
                        }]
                    }
                }
        
        # 기본 메시지 반환
        return {
            'data': [],
            'layout': {
                'title': '샘플을 선택해주세요',
                'xaxis': {'title': 'Position'},
                'yaxis': {'title': '-log10(P-value)'},
                'annotations': [{
                    'text': '샘플을 선택해 주세요',
                    'xref': 'paper', 'yref': 'paper',
                    'x': 0.5, 'y': 0.5, 'xanchor': 'center', 'yanchor': 'middle',
                    'showarrow': False, 'font': {'size': 16, 'color': 'gray'}
                }]
            }
        }
    
    '''

# =============================================================================
# 🧬 PEDIGREE BREADCRUMB CALLBACKS
# =============================================================================

@callback(
    Output('pedigree-breadcrumb', 'children'),
    [
        Input('pedigree-path-store', 'data'),
        Input('addv-apply-groups-pedi', 'data'),
        Input('fixed-vcf-item', 'data')
    ],
    prevent_initial_call=False
)
def update_pedigree_path(path_data, group_data, fixed_vcf_data):
    from dash import html

    breadcrumb = []
    path_data = path_data or []
    group_data = group_data or []
    base_name = None

    # ✅ base_name 확보 우선순위
    if fixed_vcf_data and isinstance(fixed_vcf_data, dict):
        base_name = fixed_vcf_data.get("processed_name")

    # path_store에 base entry가 있다면 그걸 우선
    base_entry = next((p for p in path_data if p.get("type") == "base"), None)
    if base_entry and not base_name:
        base_name = base_entry.get("name")

    # --- 1️⃣ base 렌더링
    if base_name:
        breadcrumb.append(html.Div([
            html.Span(base_name, className="breadcrumb-item text-dark"),
            html.Span(" (base)", style={'color': '#28a745', 'fontSize': '12px'})
        ], className="breadcrumb-item-container"))
    else:
        # ✅ 초기 값 없을 때 placeholder 표시
        return [html.Span("No pedigree loaded", style={'color': '#6c757d', 'fontStyle': 'italic'})]

    # --- 2️⃣ group-add: base 바로 뒤 표시
    if group_data:
        breadcrumb.append(html.Span(" , ", style={'margin': '0 5px', 'color': '#6c757d'}))
        for j, gname in enumerate(group_data):
            breadcrumb.append(html.Div([
                html.Span(
                    gname,
                    className="breadcrumb-item text-info",
                    style={'fontWeight': '500'}
                )
            ], className="breadcrumb-item-container"))
            if j < len(group_data) - 1:
                breadcrumb.append(html.Span(", ", style={'color': '#6c757d'}))

    # --- 3️⃣ expand
    expands = [p for p in path_data if p.get("type") == "expand" and p.get("status") != "removed"]
    for i, entry in enumerate(expands):
        name = entry["name"]
        idx = entry["index"]
        breadcrumb.append(html.Span(" > ", style={'margin': '0 5px', 'color': '#6c757d'}))
        breadcrumb.append(html.Div([
            html.Span(name, className="breadcrumb-item text-primary"),
            html.Span(f" (Expand{i+1})", style={'color': '#000', 'fontSize': '12px', 'marginLeft': '2px'}),
            html.Button(
                "×",
                id={'type': 'breadcrumb-remove', 'index': int(idx)},
                n_clicks=0,
                title=f"Remove {name}",
                style={
                    'marginLeft': '3px',
                    'border': 'none',
                    'background': 'transparent',
                    'color': '#dc3545',
                    'cursor': 'pointer',
                    'fontSize': '12px'
                }
            )
        ], className="breadcrumb-item-container"))

    return [x for x in breadcrumb if x]

'''
@callback(
    [
        Output('pedigree-path-store', 'data', allow_duplicate=True),
        Output('pedigree-breadcrumb', 'children', allow_duplicate=True),
        Output('pedigree-cytoscape', 'elements', allow_duplicate=True),
    ],
    Input({'type': 'breadcrumb-remove', 'index': ALL}, 'n_clicks'),
    [
        State('pedigree-path-store', 'data'),
        State('pedigree-path-child-store', 'data'),
        State('pedigree-cytoscape', 'elements'),
        State('addv-apply-groups-pedi', 'data'),
    ],
    prevent_initial_call=True
)
def remove_from_pedigree_path(remove_clicks, path_store, child_store, elements, group_data):
    """
    실제로 '×' 버튼이 클릭된 경우에만 expand 노드 제거.
    (새 버튼 렌더링 시 n_clicks=0 은 무시)
    """
    from dash import callback_context as ctx, no_update
    from dash.exceptions import PreventUpdate

    # --- 0️⃣ 기본 검증
    if not ctx.triggered_id or not isinstance(ctx.triggered_id, dict):
        raise PreventUpdate

    # ✅ 실제 클릭된 버튼 판별
    if not remove_clicks or all((c is None or c == 0) for c in remove_clicks):
        # 모든 n_clicks가 0이면 클릭 없음
        raise PreventUpdate

    clicked_idx = ctx.triggered_id.get("index")
    print(f"🔹 Breadcrumb remove triggered: index={clicked_idx}")

    path_store = path_store or []
    child_store = child_store or {}
    elements = elements or []
    group_data = group_data or []

    # --- 1️⃣ 대상 노드 찾기
    target_entry = next((p for p in path_store if p["index"] == clicked_idx), None)
    if not target_entry or target_entry.get("type") != "expand":
        print("⚠️ 대상이 expand 타입이 아님 — 무시")
        raise PreventUpdate

    target_name = target_entry["name"]
    print(f"🗑️ Expand 제거 요청: {target_name} (index={clicked_idx})")

    # --- 2️⃣ 상태 변경
    target_entry["status"] = "removed"

    # --- 3️⃣ 자식 노드까지 제거
    def mark_subtree_removed(parent):
        children = child_store.get(parent, [])
        for child in children:
            match = next((p for p in path_store if p["name"] == child), None)
            if match and match["status"] != "removed":
                match["status"] = "removed"
                mark_subtree_removed(child)

    mark_subtree_removed(target_name)

    # --- 4️⃣ active / base 노드 판단
    active_nodes = [p["name"] for p in path_store if p["status"] == "active"]
    base_name = next((p["name"] for p in path_store if p["type"] == "base"), None)

    if not base_name:
        print("❌ base 노드 없음 — 업데이트 불가")
        raise PreventUpdate

    if not active_nodes:
        active_nodes = [base_name]

    # --- 5️⃣ 그래프 재구성
    try:
        pedigree_app = get_pedigree_app()

        nodes, edges = pedigree_app.get_connected_nodes(base_name, 2, 2)

        for parent, childs in child_store.items():
            if parent in active_nodes:
                valid = [c for c in childs if c in active_nodes]
                for c in valid:
                    nodes.add(parent)
                    nodes.add(c)
                    edges.add((parent, c))

        new_elements = pedigree_app.create_cytoscape_elements(nodes, edges)

        # base 강조
        for el in new_elements:
            if el["data"].get("id") == base_name:
                el["classes"] = "search-result"

    except Exception as e:
        print(f"❌ 재구성 오류: {e}")
        new_elements = elements

    print(f"✅ active nodes: {active_nodes}")
    return path_store, no_update, new_elements

    '''


    

    
    # 🎯 Trait Chart Toggle 콜백 추가
@callback(
    [Output("gwas-trait-chart-collapse", "is_open"),
        Output("gwas-trait-chart-icon", "children")],
    Input("gwas-toggle-trait-chart", "n_clicks"),
    State("gwas-trait-chart-collapse", "is_open")
)
def toggle_gwas_trait_chart(n_clicks, is_open):
    """Trait Chart 접기/펼치기 토글"""
    if n_clicks:
        new_state = not is_open
        icon = "▼" if new_state else "▶"
        print(f"🎨 Trait Chart 토글: {'펼침' if new_state else '접음'}")
        return new_state, icon
    return is_open, "▼" if is_open else "▶"


# Updated 콜백: GWAS Related Filter가 적용된 경우에만 variant 섹션 표시
@callback(
    Output('gwas-variant-section', 'style'),
    Input('gwas-selected-samples-store', 'data'),#copy2gwas-selected-samples-store
    prevent_initial_call=True
)
def toggle_gwas_variant_section(selected_samples):#, selected_traits, trait_groups
    """
    샘플과 trait이 선택되면 variant 섹션(scatter plot + table)을 보여주고,
    둘 중 하나라도 없으면 숨김
    """
    print(f"🔍 toggle_gwas_variant_section 호출됨: selected_samples={selected_samples}")
    
    if selected_samples and len(selected_samples) > 0:
        # 샘플이 선택되면 2, 3번 모두 보이도록 style 변경
        show_style = {'display': 'block', 'margin-bottom': '30px'}
        print(f"   ✅ 섹션 표시: {len(selected_samples)}개 샘플 선택됨")
        return show_style#, show_style
    else:
        # 아무 샘플도 선택되지 않으면 모두 숨김
        hide_style = {'display': 'none'}
        print(f"   ❌ 섹션 숨김: 샘플 선택 없음")
        return hide_style#, hide_style
@callback(
    Output('integrated-subtrait-dropdown', 'options'),
    [
        Input('gwas-combo-store', 'data'),
        Input('analysis-tabs', 'active_tab'),
    ],
    prevent_initial_call=True
)
def update_gwas_subtrait_options(combo_json, active_tab):
    if active_tab != 'gwas-tab':
        raise PreventUpdate

    if not combo_json:
        print("[CB-subtrait] combo_json is None")
        return _default_subtrait_options()

    try:
        options = build_subtrait_options_from_combo(combo_json)
        print(f"✅ Subtrait 옵션 {len(options)}건 생성")
        return options
    except Exception as e:
        print(f"❌ Subtrait 옵션 생성 오류: {e}")
        return _default_subtrait_options()


def _default_subtrait_options():
    return [
        {'label': 'Biochemical', 'value': 'Biochemical(103)'},
        {'label': 'Biologicalprocess', 'value': 'Biologicalprocess(19)'},
        {'label': 'Plant quality', 'value': 'Plant_quality(26)'},
        {'label': 'Plantgrowth Development', 'value': 'Plantgrowth_Development(17)'},
        {'label': 'Plantmorphology', 'value': 'Plantmorphology(191)'},
        {'label': 'Plantvigor', 'value': 'Plantvigor(23)'},
        {'label': 'Stress', 'value': 'Stress(66)'},
        {'label': 'Yield', 'value': 'Yield(14)'},
    ]


@callback(
    Output('gwas-combo-store', 'data'),
    Input('gwas-selected-samples-store', 'data'), #copy2gwas-selected-samples-store
    prevent_initial_call=False
)
def update_gwas_combo_store(sample_ids):
    print(f"[CB-combo] ▶ sample_ids={sample_ids}")
    combo = get_gwas_data_for_samples_optimized(sample_ids, return_jsonable=True)
    print(f"[CB-combo] ◀ status={combo.get('status') if combo else None}, "
          f"matched={len(combo.get('matched_samples', [])) if combo else 0}, "
          f"vo={len(combo.get('variant_only', [])) if combo else 0}, "
          f"vp={len(combo.get('variant_pvalue', [])) if combo else 0}")
    return combo or {"status":"empty","requested_samples":sample_ids or [],"matched_samples":[],
                     "variant_only":[], "variant_pvalue":[], "error":"no data","combo_key":""}
    # Trait bar charts 생성
    
    
    # Trait 선택 핸들링 (차트 클릭 + span 제거 버튼)
    @callback(
        Output('gwas-selected-traits-store', 'data'),
        [Input({'type': 'gwas-trait-bar', 'subtrait': ALL}, 'clickData'),
         Input({'type': 'gwas-remove-button', 'trait': ALL}, 'n_clicks')],
        [State('gwas-selected-traits-store', 'data')]
    )
    def handle_gwas_trait_selection(click_data_list, remove_clicks, current_selected):
        ctx = callback_context
        if not ctx.triggered:
            return current_selected or []
        
        trigger = ctx.triggered[0]['prop_id']
        print(f"🎯 Trait 선택 이벤트: {trigger}")
        
        # Handle remove button clicks
        if 'gwas-remove-button' in trigger:
            trait = json.loads(trigger.split('.')[0])['trait']
            print(f"   제거: {trait}")
            return [t for t in (current_selected or []) if t != trait]
        
        # Handle chart clicks
        if click_data_list and any(click_data_list):
            for click_data in click_data_list:
                if click_data and 'points' in click_data:
                    for point in click_data['points']:
                        if 'y' in point:
                            clicked_trait = point['y']
                            current_selected = current_selected or []
                            
                            print(f"   클릭된 trait: {clicked_trait}")
                            
                            if clicked_trait in current_selected:
                                current_selected.remove(clicked_trait)
                                print(f"   제거됨: {clicked_trait}")
                            else:
                                current_selected.append(clicked_trait)
                                print(f"   추가됨: {clicked_trait}")
                            
                            return current_selected
        
        return current_selected or []

    # Selected traits span 렌더링
    @callback(
        Output('gwas-selected-traits-container', 'children'),
        Input('gwas-selected-traits-store', 'data')
    )
    def render_gwas_selected_traits(selected_traits):
        if not selected_traits:
            return []
        
        #print(f"🏷️ 선택된 traits span 생성: {selected_traits}")
        
        spans = []
        for trait in selected_traits:
            spans.append(html.Span([
                trait,
                html.Button('×', id={'type': 'gwas-remove-button', 'trait': trait}, n_clicks=0,
                            style={'margin-left': '5px', 'background': 'transparent', 'border': 'none',
                                   'color': 'red', 'cursor': 'pointer', 'font-weight': 'bold'})
            ], style={
                'display': 'inline-block', 'padding': '5px 10px', 'border': '1px solid #ccc',
                'border-radius': '4px', 'margin-right': '5px', 'margin-bottom': '5px',
                'background-color': '#e8f4fd'
            }))
        return spans
    
    # 중복 콜백 제거됨 - 이제 selected samples store만 감지
    
    # 필터/컨트롤 변경 시 플롯 업데이트 (스피닝 없음) - 통합 필터 연동
    '''
    @callback(
        Output('gwas-sample-scatter', 'figure', allow_duplicate=True),
        [Input('integrated-subtrait-dropdown', 'value'),
         Input('gwas-selected-traits-store', 'data'),
         Input('gwas-filter-states', 'data'),
         Input('gwas-input-loci-store', 'data'),
         State('gwas-pvalue-cutoff', 'value'),
         #Input('gwas-unique-position-toggle', 'value'),
         # Integrated Filter inputs
         Input('integrated-selected-traits-store', 'data'),
         Input('integrated-filter-b-state', 'data'),
         Input('integrated-filter-c-state', 'data'),
         Input('integrated-filter-d-state', 'data'),
         Input('gwas-hide-traits-checkbox', 'value')],  # hide_vals 추가
        [State('copy2gwas-selected-samples-store', 'data')],
        prevent_initial_call=True
    )
    def update_gwas_scatter_controls(selected_subtraits, selected_traits, filter_states, input_loci, 
                                   pvalue_cutoff,  integrated_traits, integrated_filter_b, 
                                   integrated_filter_c, integrated_filter_d, hide_vals, selected_samples):
        print(f"🔍 update_gwas_scatter_controls 호출됨: selected_samples={selected_samples}")
        print(f"🎛️ Scatter plot 컨트롤 업데이트 (스피닝 없음): traits={selected_traits}, integrated_traits={len(integrated_traits) if integrated_traits else 0}")
        show_unique_only=False
        # Integrated Filter의 값들을 기본 필터에 통합
        combined_traits = selected_traits or []
        if integrated_traits:
            combined_traits.extend(integrated_traits)
            combined_traits = list(set(combined_traits))  # 중복 제거
        
        # Integrated Filter B (P-value cutoff)
        if integrated_filter_b and integrated_filter_b.get('enabled', False):
            pvalue_cutoff = integrated_filter_b.get('pvalue_cutoff', pvalue_cutoff)
        
        # Integrated Filter C (MAF cutoff)
        if integrated_filter_c and integrated_filter_c.get('enabled', False):
            if not filter_states:
                filter_states = {}
            filter_states['maf'] = integrated_filter_c.get('maf_cutoff', 0.05)
        
        # Integrated Filter D (Unique position)
        if integrated_filter_d and integrated_filter_d.get('enabled', False):
            show_unique_only = True
        
        if not GWAS_AVAILABLE:
            return go.Figure().add_annotation(
                text="GWAS 모듈이 사용할 수 없습니다",
                xref="paper", yref="paper",
                x=0.5, y=0.5, xanchor='center', yanchor='middle',
                showarrow=False, font=dict(size=16, color="red")
            )
            
        if not selected_samples:
            return go.Figure().add_annotation(
                text="샘플을 선택해주세요",
                xref="paper", yref="paper",
                x=0.5, y=0.5, xanchor='center', yanchor='middle',
                showarrow=False, font=dict(size=16, color="gray")
            )
        
        try:
            gwas_result = get_gwas_data_for_samples_optimized(selected_samples)
            if not gwas_result:
                return go.Figure().add_annotation(
                    text="GWAS 데이터를 로드할 수 없습니다",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5, xanchor='center', yanchor='middle',
                    showarrow=False, font=dict(size=16, color="red")
                )
            
            gwas_module = gwas_result.get('gwas_module')
            
            # GWAS 모듈 방식 사용
            if gwas_module and hasattr(gwas_module, 'make_pvalue_scatter_for_samples'):
                maf_threshold = filter_states.get('maf') if filter_states else None
                actual_pvalue_cutoff = 10**(-pvalue_cutoff) if pvalue_cutoff else 10**(-5)
                
                fig = gwas_module.make_pvalue_scatter_for_samples(
                    selected_samples=selected_samples,
                    selected_subtraits=selected_subtraits or [],
                    maf_threshold=maf_threshold,
                    input_loci=input_loci or [],
                    pvalue_cutoff=actual_pvalue_cutoff,
                    unique_only=show_unique_only or False
                )
                
                print(f"   ✅ 컨트롤 업데이트 완료 (스피닝 없음)")
                return fig
            else:
                # 대안 방식
                fig = create_enhanced_gwas_scatter_plot(
                    gwas_result=gwas_result,
                    selected_samples=selected_samples,
                    selected_subtraits=selected_subtraits or [],
                    selected_traits=combined_traits or [],
                    filter_states=filter_states or {},
                    input_loci=input_loci or [],
                    pvalue_cutoff=pvalue_cutoff or 5,
                    show_unique_only=show_unique_only or False,
                    hide_vals=hide_vals,
                    sample_visible={sample: True for sample in selected_samples},
                    subtrait_visible={subtrait: True for subtrait in (selected_subtraits or [])}
                )
                
                print(f"   ✅ 컨트롤 업데이트 완료 (대안 방식)")
                return fig
            
        except Exception as e:
            print(f"❌ 컨트롤 업데이트 오류: {e}")
            return go.Figure().add_annotation(
                text=f"오류 발생: {str(e)}",
                xref="paper", yref="paper",
                x=0.5, y=0.5, xanchor='center', yanchor='middle',
                showarrow=False, font=dict(size=14, color="red")
            )
    '''
    '''
    # Download container visibility
    @callback(
        Output('gwas-download-container', 'style'),
        Input('copy2gwas-selected-samples-store', 'data'),
        prevent_initial_call=True  # 초기 로딩도 허용
    )
    def control_gwas_download_visibility(selected_samples):
        if selected_samples:
            return {'display': 'block'}
        else:
            return {'display': 'none'}
    
    # Locus table 업데이트 (기본적으로 표시되도록 수정)
    @callback(
        [Output('gwas-locus-table', 'columns'),
         Output('gwas-locus-table', 'data')],
        [Input('copy2gwas-selected-samples-store', 'data'),  # 샘플만으로도 기본 테이블 표시
         Input('gwas-selected-traits-store', 'data'),
         Input('integrated-subtrait-dropdown', 'value'),
         Input('gwas-filter-states', 'data'),
         Input('gwas-input-loci-store', 'data'),
         Input('integrated-trait-groups-store', 'data'),
         Input('integrated-secondary-filter-enabled', 'data')],
        prevent_initial_call=True  # 초기 로딩 방지
    )
    def update_gwas_locus_table(selected_samples, selected_traits, selected_subtraits, filter_states, input_loci, trait_groups, secondary_enabled):
        print(f"📋 Locus table 업데이트: samples={selected_samples}, traits={selected_traits}")
        
        if not GWAS_AVAILABLE:
            print(f"   ❌ Missing requirements: GWAS_AVAILABLE={GWAS_AVAILABLE}")
            return [], []
        
        if not selected_samples:
            print(f"   ⏳ No samples selected yet - showing empty table with headers")
            # Return empty table with basic column structure
            empty_columns = [
                {'name': 'Chromosome', 'id': 'Chromosome', 'type': 'text'},
                {'name': 'Position', 'id': 'Position', 'type': 'numeric'},
                {'name': 'P-value', 'id': 'P-value', 'type': 'numeric'},
                {'name': 'Trait', 'id': 'Trait', 'type': 'text'},
                {'name': 'Subtrait', 'id': 'Subtrait', 'type': 'text'},
            ]
            return empty_columns, []
        
        try:
            gwas_result = get_gwas_data_for_samples_optimized(selected_samples)
            if not gwas_result:
                return [], []
            
            gwas_module = gwas_result['gwas_module']
            # GWASModule의 update_locus_table_data 메서드 사용
            columns, data = gwas_module.update_locus_table_data(
                selected_traits or [], selected_samples, selected_subtraits, filter_states or {}, input_loci or []
            )
            
            # 2차 필터링이 활성화된 경우 그룹 컬럼 추가
            if secondary_enabled and trait_groups and trait_groups.get('groups'):
                # 그룹 정보 파싱
                groups = trait_groups['groups']
                group_trait_map = {}  # trait -> group_number 매핑
                
                for group_name, traits_in_group in groups.items():
                    if traits_in_group is not None and len(traits_in_group) > 0:
                        group_number = group_name.replace('Group ', '').replace('group', '').strip()
                        for trait in traits_in_group:
                            group_trait_map[trait] = group_number
                
                # 그룹 컬럼 추가
                if group_trait_map:
                    # GROUP 컬럼 추가 (단일 컬럼)
                    enhanced_columns = columns.copy()
                    enhanced_columns.append({
                        'name': 'GROUP',
                        'id': 'GROUP',
                        'type': 'text'
                    })
                    
                    # 데이터에 그룹 정보 추가
                    enhanced_data = []
                    for row in data:
                        enhanced_row = row.copy()
                        trait = row.get('Trait', '')
                        
                        # 트레이트가 어떤 그룹에 속하는지 확인
                        group_num = group_trait_map.get(trait, '')
                        enhanced_row['GROUP'] = group_num
                        enhanced_data.append(enhanced_row)
                    
                    print(f"   ✅ GROUP 컬럼 추가: {len(enhanced_columns)}열, {len(enhanced_data)}행")
                    return enhanced_columns, enhanced_data
            
            print(f"   ✅ 테이블 생성: {len(columns)}열, {len(data)}행")
            return columns, data
            
        except Exception as e:
            print(f"❌ 테이블 업데이트 오류: {e}")
            import traceback
            traceback.print_exc()
            return [], []
    
    # Download all loci
    @callback(
        Output("gwas-download-all-loci", "data"),
        Input("gwas-download-all-btn", "n_clicks"),
        [State('copy2gwas-selected-samples-store', 'data'),
         State('integrated-subtrait-dropdown', 'value'),
         State('gwas-filter-states', 'data'),
         State('gwas-input-loci-store', 'data')],
        prevent_initial_call=True
    )
    def download_gwas_all_loci(n_clicks, selected_samples, selected_subtraits, filter_states, input_loci):
        print(f"💾 Download all loci: {n_clicks}")
        if n_clicks and selected_samples and GWAS_AVAILABLE:
            try:
                gwas_result = get_gwas_data_for_samples_optimized(selected_samples)
                if gwas_result:
                    gwas_module = gwas_result['gwas_module']
                    return gwas_module._handle_download_all_loci(selected_samples, selected_subtraits, filter_states or {}, input_loci or [])
            except Exception as e:
                print(f"❌ Download 오류: {e}")
        return no_update
'''

# 🚩 Reset View button callback (backup method)
@callback(
    [
        Output('pedigree-cytoscape', 'elements', allow_duplicate=True),
        Output('pedigree-cytoscape', 'layout', allow_duplicate=True),
        Output('reset-view-trigger', 'data'),
        Output('pedigree-path-store', 'data', allow_duplicate=True),
        Output('pedigree-path-child-store', 'data', allow_duplicate=True),
        Output('selected-nodes-store', 'data', allow_duplicate=True),
        Output('pedigree-cytoscape', 'stylesheet', allow_duplicate=True),
    ],
    Input('btn-reset', 'n_clicks'),
    [
        State('fixed-vcf-item', 'data'),
        State('fixed-vcf-item-nopedi', 'data'),  # 🆕 추가
        State('addv-apply-groups-pedi', 'data'),
        State('addv-apply-groups-nopedi', 'data'),
    ],
    prevent_initial_call=True
)
def handle_reset_view_full(n_clicks, fixed_vcf,fixed_vcf_item_nopedi, group_pedi, group_nopedi):
    """Reset View: 초기 pedigree + nopedi 구성요소 재생성"""
    from dash.exceptions import PreventUpdate
    import math

    if not n_clicks:
        raise PreventUpdate

    print(f"🔄 Full Reset View triggered ({n_clicks})")

        # ✅ base node 병합
    fixed_nodes = []

    if fixed_vcf and isinstance(fixed_vcf, dict):
        fixed_nodes.append({
            "id": str(fixed_vcf.get("processed_name") or ""),
            "vcf_id": str(fixed_vcf.get("status") or ""),
        })

    if fixed_vcf_item_nopedi and isinstance(fixed_vcf_item_nopedi, dict):
        fixed_nodes.append({
            "id": str(fixed_vcf_item_nopedi.get("variety_id") or ""),
            "vcf_id": str(fixed_vcf_item_nopedi.get("status") or ""),
        })

    # ✅ basenode 결정 (둘 중 하나만 사용됨)
    basenode = None
    if fixed_nodes:
        basenode = next((n.get("id") for n in fixed_nodes if n.get("id")), None)



    pedigree_app = get_pedigree_app()
    if pedigree_app is None:
        print("❌ get_pedigree_app() returned None. Abort reset.")
        raise PreventUpdate

    elements = []
    try:
        # -----------------------------
        # ① Base pedigree
        # -----------------------------
        fixed_name = None
        if fixed_vcf and isinstance(fixed_vcf, dict):
            fixed_name = fixed_vcf.get("processed_name") or fixed_vcf.get("variety_id")
        if fixed_name:
            base_nodes, base_edges = pedigree_app.get_connected_nodes(fixed_name, 2, 2)
            base_elements = pedigree_app.create_cytoscape_elements(base_nodes, base_edges)
            elements.extend(base_elements)
            print(f"🟩 Base pedigree built for {fixed_name} ({len(base_nodes)} nodes)")

        elif fixed_vcf_item_nopedi and isinstance(fixed_vcf_item_nopedi, dict):
            print("🌾 NOPEDI case detected — generating standalone node")

            node_id = str(fixed_vcf_item_nopedi.get("variety_id") or "").strip()
            vcf_status = str(fixed_vcf_item_nopedi.get("status") or "No VCF").strip()

            has_vcf = False if vcf_status in ("", "No VCF", None) else True
            color_class = get_color_class(has_vcf, True)
            elements.append({
                "data": {
                    "id": node_id,
                    "label": node_id,
                    "type": "nopedi",
                    "it_number": node_id,
                    "vcf_id": vcf_status,
                    "vcf_status": vcf_status,
                    "has_vcf": has_vcf,
                    "has_it": True,
                    "color_class": color_class,
                },
                "position": {"x": 0, "y": 0},
                "classes": "nopedi-node",
                "style": {
                    "shape": "rectangle",
                    "width": 90,
                    "height": 30,
                    "font-size": "11px",
                    "text-valign": "center",
                    "text-halign": "center",
                    "border-width": 1.5,
                    # ✅ 시각 효과
                    "shadow-color": "#7f8c8d",
                    "shadow-blur": 10,
                    "shadow-opacity": 0.6,
                    "shadow-offset-x": 2,
                    "shadow-offset-y": 2,
                    "underlay-color": "#bdc3c7",
                    "underlay-opacity": 0.8,
                    "underlay-padding": 3,
                },
            })
                
        # -----------------------------
        # ② Group_pedi 추가
        # -----------------------------
        if group_pedi:
            for g in group_pedi:
                gn, ge = pedigree_app.get_connected_nodes(g, 2, 2)
                g_elem = pedigree_app.create_cytoscape_elements(gn, ge)
                elements.extend(g_elem)
                print(f"🟦 Added group lineage for {g} ({len(gn)} nodes)")

        # -----------------------------
        # ③ NOPEDI 반영
        # -----------------------------
        if isinstance(group_nopedi, dict) and "records" in group_nopedi:
            nopedi_records = group_nopedi.get("records", [])
            total = len(nopedi_records)
            print(f"🟧 Adding {total} NOPEDI nodes (hourglass layout)")

            merged_elements = elements.copy()
            existing_node_ids = {e["data"]["id"] for e in merged_elements if "data" in e and "id" in e["data"]}

            def _safe_key(rec):
                name = rec.get("name")
                it = rec.get("it_number")
                name = name if isinstance(name, str) else ("" if name is None else str(name))
                it = it if isinstance(it, str) else ("" if it is None else str(it))
                name, it = name.strip(), it.strip()
                if name.lower() in ("none", "null", "nan", ""):
                    name = ""
                if it.lower() in ("none", "null", "nan", ""):
                    it = ""
                return name or it

            batch_size = 4
            patterns, pending_upper = [], []
            for i in range(0, total, batch_size):
                subset = nopedi_records[i:i + batch_size]
                if not pending_upper:
                    pending_upper = subset
                    patterns.append(("side", subset))
                else:
                    patterns.append(("hourglass", (pending_upper, subset)))
                    pending_upper = []

            base_x, base_y = 0, 1000
            x_spacing, y_spacing = 100, 150
            side_offset, anchor_gap = 250, 300
            anchor_index, anchor_y = 1, base_y
            existing_anchor = None

            for ptype, content in patterns:
                if ptype == "side":
                    subset = content
                    anchor_id = f"nopedi_anchor_{anchor_index}"
                    anchor_index += 1

                    merged_elements.append({
                        "data": {"id": anchor_id, "label": f"NOPEDI SIDE-{anchor_index}"},
                        "classes": "nopedi-anchor",
                        "position": {"x": base_x + side_offset, "y": anchor_y},
                        "style": {
                            "background-color": "#ffffff",
                            "border-color": "#555",
                            "border-style": "dashed",
                            "border-width": 2,
                            "shape": "ellipse",
                            "width": 80,
                            "height": 80,
                            "opacity": 0,
                            #"label": f"SIDE-{anchor_index}",
                            "font-size": "10px",
                        },
                    })

                    for i, rec in enumerate(subset):
                        node_key = _safe_key(rec)
                        it_number = rec.get("it_number", "")
                        has_vcf = rec.get("has_vcf", False)
                        has_pheno = rec.get("has_pheno", False)
                        vcf_id = rec.get("vcf_id", "No VCF")
                        color_class = get_color_class(has_vcf, has_pheno)
                        if not node_key or node_key in existing_node_ids:
                            continue
                        x = base_x + side_offset + (i - len(subset)/2) * x_spacing
                        y = anchor_y
                        merged_elements.append({
                            "data": {"id": node_key, "label": node_key, "type": "nopedi",
                             "it_number": it_number,'vcf_id': vcf_id,'vcf_status': vcf_id ,
                             "has_vcf": has_vcf, "has_it": has_pheno, "color_class": color_class},
                            "position": {"x": x, "y": y},
                            "classes": "nopedi-base-node",
                            "style": {
                                #"background-color": "#f39c12",
                                "shape": "rectangle",
                                "font-size": "10px",
                                "text-valign": "center",
                                "text-halign": "center",
                            },
                        })
                        merged_elements.append({
                            "data": {"source": anchor_id, "target": node_key},
                            "classes": "nopedi-link",
                            "style": {"line-style": "dotted", "width": 1, "line-color": "#ccc", "opacity": 0},
                        })
                        existing_node_ids.add(node_key)
                    anchor_y += anchor_gap

                elif ptype == "hourglass":
                    upper, lower = content
                    anchor_id = f"nopedi_anchor_{anchor_index}"
                    anchor_index += 1
                    merged_elements.append({
                        "data": {"id": anchor_id, "label": f"NOPEDI HG-{anchor_index}"},
                        "classes": "nopedi-anchor",
                        "position": {"x": base_x, "y": anchor_y},
                        "style": {
                            "background-color": "#ffffff",
                            "border-color": "#555",
                            "border-style": "dashed",
                            "border-width": 2,
                            "shape": "ellipse",
                            "width": 90,
                            "height": 90,
                            "opacity": 0,
                            "label": f"HG-{anchor_index}",
                            "font-size": "10px",
                        },
                    })
                    for i, rec in enumerate(upper):
                        node_key = _safe_key(rec)
                        it_number = rec.get("it_number", "")
                        has_vcf = rec.get("has_vcf", False)
                        has_pheno = rec.get("has_pheno", False)
                        vcf_id = rec.get("vcf_id", "No VCF")
                        color_class = get_color_class(has_vcf, has_pheno)
                        if not node_key or node_key in existing_node_ids:
                            continue
                        x = base_x + (i - len(upper)/2) * x_spacing
                        y = anchor_y - y_spacing
                        merged_elements.append({
                            "data": {"id": node_key, "label": node_key, "type": "nopedi", 
                            "it_number": it_number,'vcf_id': vcf_id,'vcf_status': vcf_id ,
                            "has_vcf": has_vcf, "has_it": has_pheno, "color_class": color_class},
                            "position": {"x": x, "y": y},
                            "classes": "nopedi-node",
                            "style": {
                                #"background-color": "#f39c12",
                                "shape": "rectangle",
                                "font-size": "10px",
                            },
                        })
                        merged_elements.append({
                            "data": {"source": node_key, "target": anchor_id},
                            "classes": "nopedi-link",
                            "style": {"line-style": "dotted", "width": 1, "line-color": "#ccc", "opacity": 0},
                        })
                        existing_node_ids.add(node_key)
                    for i, rec in enumerate(lower):
                        node_key = _safe_key(rec)
                        it_number = rec.get("it_number", "")
                        has_vcf = rec.get("has_vcf", False)
                        has_pheno = rec.get("has_pheno", False)
                        vcf_id = rec.get("vcf_id", "No VCF")
                        color_class = get_color_class(has_vcf, has_pheno)
                        if not node_key or node_key in existing_node_ids:
                            continue
                        x = base_x + (i - len(lower)/2) * x_spacing
                        y = anchor_y + y_spacing
                        merged_elements.append({
                            "data": {"id": node_key, "label": node_key, "type": "nopedi", 
                            "it_number": it_number,'vcf_id': vcf_id,'vcf_status': vcf_id ,
                            "has_vcf": has_vcf, "has_it": has_pheno, "color_class": color_class},
                            "position": {"x": x, "y": y},
                            "classes": "nopedi-node",
                            "style": {
                                #"background-color": "#f39c12",
                                "shape": "rectangle",
                                "font-size": "10px",
                            },
                        })
                        merged_elements.append({
                            "data": {"source": anchor_id, "target": node_key},
                            "classes": "nopedi-link",
                            "style": {"line-style": "dotted", "width": 1, "line-color": "#ccc", "opacity": 0},
                        })
                        existing_node_ids.add(node_key)
                    if existing_anchor:
                        merged_elements.append({
                            "data": {"source": existing_anchor, "target": anchor_id},
                            "classes": "anchor-bridge",
                            "style": {"line-style": "dashed", "width": 2, "line-color": "#999", "opacity": 0},
                        })
                    existing_anchor = anchor_id
                    anchor_y += anchor_gap

            elements = merged_elements
            print(f"🟧 NOPEDI merged ({len(elements)} total elements)")

        # -----------------------------
        # ④ 중복 제거
        # -----------------------------
        unique_ids, deduped = set(), []
        for e in elements:
            key = (e["data"].get("source"), e["data"].get("target")) if "source" in e["data"] else e["data"]["id"]
            if key not in unique_ids:
                deduped.append(e)
                unique_ids.add(key)
        elements = deduped

    except Exception as e:
        print(f"❌ Reset build error: {e}")
        raise PreventUpdate

    reset_layout = {
        'name': 'dagre',
        'rankDir': 'TB',
        'ranker': 'network-simplex',
        'animate': True,
        'animationDuration': 800,
        'fit': True,
        'padding': 60,
        'spacingFactor': 1.1,
        'nodeDimensionsIncludeLabels': True,
        'rankSep': 100,
        'nodeSep': 80,
        'edgeSep': 40,
        'avoidOverlap': True,
        'refresh': n_clicks
    }
    style = get_default_stylesheet2(basenode)

    print(f"✅ Pedigree fully reset ({len(elements)} elements incl. NOPEDI)")
    return elements, reset_layout, n_clicks, [], {}, [], style


'''
@callback(
    [Output('pedigree-cytoscape', 'elements', allow_duplicate=True),
     Output('selected-node-store', 'data', allow_duplicate=True),
     Output('multi-selected-nodes', 'data', allow_duplicate=True)],
    Input('search-button', 'n_clicks'),
    [State('search-input', 'value'),
     State('multi-selected-nodes', 'data')],
    prevent_initial_call=True
)
def handle_search(n_clicks, search_value, current_selected):
    """통합된 검색 콜백 (중복 제거) - 다중 선택 목록에 자동 추가"""
    if not n_clicks or not search_value or not GWAS_AVAILABLE:
        return [], None, current_selected or []
    
    try:
        print(f"🔍 검색 실행: {search_value}")
        pedigree_app = get_pedigree_app()
        nodes, edges = pedigree_app.get_connected_nodes(search_value, 2, 2)
        elements = pedigree_app.create_cytoscape_elements(nodes, edges)
        
        # 검색된 노드 하이라이트
        for element in elements:
            if ('data' in element and 'id' in element['data'] and 
                element['data']['id'] == str(search_value)):
                element['classes'] = 'search-result'
        
        # 검색된 노드를 다중 선택 목록에 추가 (중복 방지)
        updated_selected = list(current_selected or [])
        if search_value not in updated_selected:
            updated_selected.append(search_value)
            print(f"   🎯 검색된 노드 '{search_value}' 다중 선택 목록에 추가")
        
        print(f"   생성된 elements: {len(elements)}개")
        print(f"   업데이트된 선택 목록: {updated_selected}")
        print(f"   🔍 검색 완료 - VCF 자동 추출 및 Fixed sample 자동 체크 예정")
        return elements, search_value, updated_selected
        
    except Exception as e:
        print(f"❌ Error in search: {e}")
        return [], None, current_selected or []
'''

'''
@callback(
    [Output('tooltip', 'children'),
     Output('tooltip', 'style',allow_duplicate=True)],
    Input('pedigree-cytoscape', 'mouseoverNodeData'),
    prevent_initial_call=True
)
def update_tooltip(over_data):
    #print(over_data)
    """🎯 Dash callback 방식 tooltip - mouseover 만 처리"""
    
    # mouseover인 경우에만 툴팁 표시
    if over_data:
        print(f"🔍 Tooltip mouseover: {over_data.get('id', 'Unknown')}")
        
        # 노드 데이터에서 정보 추출
        variety_id = over_data.get('id', 'Unknown')
        variety_label = over_data.get('label', variety_id)
        
        # 노드 위치 정보 가져오기 (position이 있다면 사용)
        position = over_data.get('position', {})
        
        # 노드 데이터에서 정보 추출
        parents_data = over_data.get('parents', [])
        children_data = over_data.get('children', [])
        it_number = over_data.get('it_number', 'No information')
        vcf_status = over_data.get('vcf_status', 'No VCF')
        has_it = over_data.get('has_it', False)
        has_vcf = over_data.get('has_vcf', False)
        
        # Parents 정보 처리
        if parents_data and len(parents_data) > 0:
            parents_display = parents_data if isinstance(parents_data, list) else [parents_data]
            if len(parents_display) > 2:
                parents_info = f"{len(parents_display)} parents: {', '.join(parents_display[:2])}..."
            else:
                parents_info = f"{len(parents_display)} parent(s): {', '.join(parents_display)}"
        else:
            parents_info = "No parents"
        
        # Children 정보 처리  
        if children_data and len(children_data) > 0:
            children_display = children_data if isinstance(children_data, list) else [children_data]
            if len(children_display) > 2:
                children_info = f"{len(children_display)} children: {', '.join(children_display[:2])}..."
            else:
                children_info = f"{len(children_display)} children: {', '.join(children_display)}"
        else:
            children_info = "No children"
        
        # IT 정보 처리
        if has_it and it_number != "No information":
            phenotype_info = f"✅ IT: {it_number}"
        else:
            phenotype_info = "❌ No IT number"
        
        # VCF 정보 처리
        if has_vcf and vcf_status != "No VCF":
            vcf_info = f"✅ VCF: {vcf_status}"
        else:
            vcf_info = "❌ No VCF data"
        
        # Tooltip 내용 생성
        tooltip_content = [
            html.Div([
                # Header with close button
                html.Div([
                    html.Div([
                        html.Strong(variety_label, style={'color': '#fff', 'fontSize': '16px', 'marginBottom': '8px', 'display': 'block'}),
                        html.Small(f"ID: {variety_id}", style={'color': '#ccc', 'fontSize': '11px'})
                    ], style={'flex': '1'}),
                    html.Button(
                        "×",
                        id={'type': 'tooltip-remove-btn', 'index': variety_id},
                        className="btn-close btn-sm",
                        title=f"Remove {variety_id} from selection",
                        style={
                            'background': 'none',
                            'border': 'none',
                            'color': '#fff',
                            'fontSize': '18px',
                            'padding': '0',
                            'marginLeft': '8px',
                            'cursor': 'pointer',
                            'opacity': '0.7',
                            'lineHeight': '1'
                        }
                    )
                ], style={
                    'display': 'flex',
                    'alignItems': 'flex-start',
                    'borderBottom': '1px solid #555', 
                    'paddingBottom': '8px', 
                    'marginBottom': '8px'
                }),
                
                # Pedigree info
                html.Div([
                    html.Div([
                        html.I(className="fas fa-level-up-alt", style={'marginRight': '6px', 'color': '#ffc107'}),
                        html.Span(parents_info, style={'fontSize': '11px'})
                    ], style={'marginBottom': '4px'}),
                    html.Div([
                        html.I(className="fas fa-level-down-alt", style={'marginRight': '6px', 'color': '#17a2b8'}),
                        html.Span(children_info, style={'fontSize': '11px'})
                    ], style={'marginBottom': '4px'})
                ], style={'marginBottom': '8px'}),
                
                # Data availability
                html.Div([
                    html.Div([
                        html.I(className="fas fa-chart-bar", style={'marginRight': '6px', 'color': '#28a745'}),
                        html.Span(phenotype_info, style={'fontSize': '11px'})
                    ], style={'marginBottom': '4px'}),
                    html.Div([
                        html.I(className="fas fa-dna", style={'marginRight': '6px', 'color': '#6f42c1'}),
                        html.Span(vcf_info, style={'fontSize': '11px'})
                    ], style={'marginBottom': '6px'})
                ]),
                
                # Footer
                html.Div([
                    html.I(className="fas fa-mouse-pointer", style={'marginRight': '4px', 'color': '#6c757d'}),
                    html.Span("Click to expand/select", style={'fontSize': '9px', 'fontStyle': 'italic', 'color': '#aaa'})
                ])
            ])
        ]
        
        # Tooltip 스타일 - 노드를 따라 움직이는 상대적 위치
        tooltip_style = {
            'position': 'absolute',  # fixed에서 absolute로 변경
            'backgroundColor': 'rgba(0, 0, 0, 0.92)',
            'color': 'white',
            'padding': '12px',
            'borderRadius': '8px',
            'fontSize': '12px',
            'display': 'block',
            'visibility': 'visible',
            'opacity': '1',
            'zIndex': 10000,
            'pointerEvents': 'auto',  # tooltip의 버튼 클릭을 허용하도록 변경
            'maxWidth': '280px',
            'minWidth': '220px',
            'border': '1px solid #28a745',
            'boxShadow': '0 4px 12px rgba(0,0,0,0.4)',
            'lineHeight': '1.4',
            'wordWrap': 'break-word',
            'bottom': '300px',        # ✅ Cytoscape 아래에서 20px 위
            'left': '20px',   
        }

        
        return tooltip_content, tooltip_style
    
    # 기본값: 툴팁 숨김
    return "", {'display': 'none', 'visibility': 'hidden', 'opacity': '0'}
'''

@callback(
    [Output('tooltip', 'children', allow_duplicate=True),
     Output('tooltip', 'style', allow_duplicate=True)],
    Input({'type': 'tooltip-remove-btn', 'index': ALL}, 'n_clicks'),
    prevent_initial_call=True
)
def close_tooltip(remove_clicks):
    """Tooltip의 x버튼 클릭 시 tooltip만 닫기 (노드는 제거하지 않음)"""
    try:
        from dash import ctx, no_update
        
        # 클릭이 없으면 변경 없음
        if not any(remove_clicks) or not ctx.triggered:
            return no_update, no_update
        
        print(f"🔍 Tooltip close clicked")
        
        # Tooltip 숨기기
        tooltip_style = {'display': 'none', 'visibility': 'hidden', 'opacity': '0'}
        
        return "", tooltip_style
    
    except Exception as e:
        print(f"❌ Tooltip close 오류: {e}")
        from dash import no_update
        return no_update, no_update



# Helper function to validate which nodes have phenotype data
# def validate_nodes_with_data(node_ids):
# 함수 제거됨 - available-nodes-store가 이미 유효한 노드들만 제공하므로 불필요

# Pedi 옵션 별도 콜백 제거됨 - 메인 콜백에서 조건부 처리로 변경

# 탭 콘텐츠 업데이트 (확장된 노드 변경에도 반응)
'''
@callback(
    Output('tab-content', 'children'),
    [
        Input('analysis-tabs', 'active_tab'),
        Input('applied-varieties-store', 'data'),     # Apply된 선택 품종들
        Input('phenotype-data-store', 'data'),        # Apply로 전달된 표현형 샘플
        Input('child-nodes-store', 'data'),           # 확장으로 추가된 노드들
        Input('pedigree-elements-store', 'data'),     # (사용은 안 하지만 트리거 용도로 유지)
        Input('search-selected-variety-name', 'data') # 검색 선택 품종 (트리거 용)
    ],
    [
        State('fixed-vcf-item', 'data'),              # regular 초기 품종 (검색품종)
        State('fixed-vcf-item-nopedi', 'data')        # nopedi 초기 품종 (IT 번호)
    ],
    prevent_initial_call=True  # 초기 렌더에서 바로 실행되도록 권장
)
def update_tab_content(active_tab,
                       applied_nodes,
                       phenotype_data,
                       child_nodes,
                       pedigree_elements,
                       search_selected_variety,
                       fixed_vcf_value,
                       fixed_vcf_nopedi):
    from dash import ctx, no_update
    trig = ctx.triggered[0]['prop_id'] if ctx.triggered else ''

    # ✅ 초기 렌더 직후: active_tab 트리거만 들어왔고, 다른 입력이 모두 비었으면 건드리지 않음
    if trig == 'analysis-tabs.active_tab' and active_tab == 'phenotype-tab':
        if not any([applied_nodes, phenotype_data, child_nodes, pedigree_elements,
                    search_selected_variety, fixed_vcf_value, fixed_vcf_nopedi]):
            return no_update
    print('트리거됬나요?')

    # ---------------- 공통 상태 파싱 ----------------
    applied_data = applied_nodes or []
    #print(f"applied_data: {applied_data}")
    selected_nodes = [
        d.get('it_number') for d in applied_data
        if isinstance(d, dict) and d.get('it_number')
    ]
    #print(f"selected_nodes: {selected_nodes}")
    phenotype_samples = phenotype_data or []
    #print(f"phenotype_samples: {phenotype_samples}")
    #print(f"pedigree_elements: {pedigree_elements}")

    is_nopedi_case   = bool(isinstance(fixed_vcf_nopedi, dict) and fixed_vcf_nopedi.get('variety_id'))
    is_regular_case  = bool(isinstance(fixed_vcf_value,   dict) and fixed_vcf_value.get('processed_name'))
    #print(f"is_nopedi_case: {is_nopedi_case}")
    #print(f"is_regular_case: {is_regular_case}")

    # ---------------- PHENOTYPE 탭 ----------------
    if active_tab == "phenotype-tab":
        # 분석 대상 노드 합성 (nopedi 기본 + regular 기본 + applied + child + phenotype-data)
        analysis_nodes = []

        # regular 초기 품종: 가능하면 IT 번호(variety_id) 우선, 없으면 processed_name
        if isinstance(fixed_vcf_value, dict):
            base_reg_it = fixed_vcf_value.get('variety_id')
            base_reg_nm = fixed_vcf_value.get('processed_name')
            base_regular = base_reg_it or base_reg_nm
            if base_regular:
                analysis_nodes.append(str(base_regular))

        # nopedi 초기 품종: IT 번호
        if isinstance(fixed_vcf_nopedi, dict):
            base_nopedi = fixed_vcf_nopedi.get('variety_id')
            if base_nopedi:
                analysis_nodes.append(str(base_nopedi))

        # applied 선택 노드들
        if selected_nodes:
            analysis_nodes.extend([str(x) for x in selected_nodes])

        # child 노드들
        if child_nodes:
            analysis_nodes.extend([str(x) for x in child_nodes])

        # phenotype-data-store 샘플들
        for s in phenotype_samples:
            if isinstance(s, dict):
                sid = s.get('id') or s.get('variety_id') or s.get('name')
                if sid:
                    analysis_nodes.append(str(sid))

        # 중복 제거 + 빈 문자열 제거
        analysis_nodes = [x for x in analysis_nodes if x]
        analysis_nodes = list(dict.fromkeys(analysis_nodes))
        print(f"analysis_nodes: {analysis_nodes}")
        # PHENO_DF 확인
        try:
            _empty = (PHENO_DF is None or PHENO_DF.empty)
        except NameError:
            _empty = True

        if _empty:
            return html.Div([
                dbc.Alert([
                    html.H5("❌ No Phenotype Data", className="alert-heading"),
                    html.P("PHENO_DF is empty or not loaded. Ensure `phenotype_data_en.csv` is loaded to PHENO_DF.")
                ], color="danger")
            ])

        # 분석 대상이 아직 없으면 안내 메시지
        if not analysis_nodes:
            return html.Div([
                dbc.Alert([
                    html.H5("🌾 Analysis System Ready", className="alert-heading"),
                    html.P("Select varieties in selection mode and click 'Additional Analysis' button.")
                ], color="info")
            ])

        # 통합 패널 생성 (단일 데이터 소스: PHENO_DF)
        return create_phenotype_panel_unified(
            selected_nodes=analysis_nodes,
            pheno_df=PHENO_DF,
            table_enabled=True
        )

    # ---------------- GWAS 탭 ----------------
    elif active_tab == "gwas-tab":
        vcf_status = None
        if is_nopedi_case and isinstance(fixed_vcf_nopedi, dict):
            vcf_status = fixed_vcf_nopedi.get('status')
        elif isinstance(fixed_vcf_value, dict):
            vcf_status = fixed_vcf_value.get('status')

        if vcf_status and vcf_status != 'No VCF':
            return create_gwas_analysis_content(vcf_status)

        msg   = "No VCF data found for this nopedi variety." if is_nopedi_case else "VCF data is required for GWAS analysis."
        color = "warning" if is_nopedi_case else "danger"
        return html.Div([
            dbc.Alert([
                html.H5("❌ No VCF data", className="alert-heading"),
                html.P(msg)
            ], color=color)
        ])

    # ---------------- 그 외 탭 ----------------
    return no_update
    '''

# GWAS 탭 내용 업데이트 (VCF 값 변경 시에만)
'''
@callback(
    Output('tab-content', 'children', allow_duplicate=True),
    [Input('current-vcf-values', 'data')],
    [State('analysis-tabs', 'active_tab'),
     State('multi-selected-nodes', 'data')],
    prevent_initial_call=True
)
def update_gwas_tab_on_vcf_change(vcf_values, active_tab, selected_nodes):
    """VCF 값이 변경되었을 때만 GWAS 탭 내용 업데이트 - 제한적으로 실행"""
    # 체크박스만으로 관리하기 위해 이 콜백을 거의 비활성화
    print(f"🔄 VCF 탭 업데이트 시도 (제한모드): vcf_values={vcf_values}, active_tab={active_tab}")
    
    # 아주 특별한 경우에만 탭 내용을 업데이트 (예: 완전히 새로운 샘플 추가)
    # 일반적인 체크박스 변경은 이 콜백을 트리거하지 않도록 함
    raise PreventUpdate
'''


setup_gwas_callbacks()

# ==============================================
# 📋 INTEGRATED FILTER CALLBACKS (통합 필터 시스템)
# ==============================================

def setup_integrated_filter_callbacks():
    """통합 필터 시스템 콜백 설정"""
    from dash.exceptions import PreventUpdate
    import dash
    
    # 통합 필터 섹션 표시 여부 제어
    '''
    @callback(
        Output('integrated-filter-section', 'style'),
        Input('copy2gwas-selected-samples-store', 'data')
    )
    def control_integrated_filter_visibility(selected_samples):
        if selected_samples:
            return {'display': 'block', 'margin-bottom': '15px', 'max-width': '100%'}
        return {'display': 'none', 'margin-bottom': '15px', 'max-width': '100%'}
    '''
    '''
    # Integrated Filter D (Unique) - 샘플 개수 제한 콜백 (warning + disabled 처리)
    @callback(
        [Output('integrated-filter-d-sample-warning', 'children'),
         Output('integrated-unique-position-enabled', 'disabled')],
        [Input('copy2gwas-selected-samples-store', 'data')]
    )
    def update_filter_d_sample_restriction(selected_samples):
        """Filter D (Unique)에 대한 샘플 개수 제한 적용"""
        if not selected_samples or len(selected_samples) < 2:
            # 샘플이 2개 미만인 경우
            warning_content = dbc.Alert([
                html.I(className="fas fa-exclamation-triangle", style={'margin-right': '8px'}),
                f"Warning: Need at least 2 samples for unique filtering. Currently selected: {len(selected_samples) if selected_samples else 0}"
            ], color="warning", style={'font-size': '0.9em', 'margin-bottom': '10px'})
            
            return warning_content, True  # disabled=True
        else:
            # 샘플이 2개 이상인 경우
            info_content = dbc.Alert([
                html.I(className="fas fa-check-circle", style={'margin-right': '8px'}),
                f"Ready: {len(selected_samples)} samples selected - unique filtering available"
            ], color="success", style={'font-size': '0.9em', 'margin-bottom': '10px'})
            
            return info_content, False  # disabled=False
    '''
    '''
    # 탭 가시성 제어 콜백
    @callback(
        [Output('integrated-filter-tab-a', 'style'),
         Output('integrated-filter-tab-b', 'style'),
         Output('integrated-filter-tab-c', 'style'),
         Output('integrated-filter-tab-d', 'style')],
        Input('filter-tabs-visibility', 'data'),
        prevent_initial_call=True
    )
    def control_filter_tabs_visibility(visibility_data):
        """Filter A 상태에 따라 탭 가시성 제어"""
        if not visibility_data:
            # 기본 상태: 모든 탭 보이기
            default_style = {'width': '170px', 'height': '60px', 'margin': '3px', 'display': 'inline-block'}
            return [default_style] * 4
        
        styles = []
        default_style = {'width': '170px', 'height': '60px', 'margin': '3px', 'display': 'inline-block'}
        hidden_style = {'width': '170px', 'height': '60px', 'margin': '3px', 'display': 'none'}
        
        for tab in ['a', 'b', 'c', 'd']:
            if visibility_data.get(tab, True):
                styles.append(default_style)
            else:
                styles.append(hidden_style)
        
        return styles
        '''
    
    # Enhanced Filter Summary 업데이트 콜백

     
     # 필터 탭 상태 표시 업데이트 콜백
    '''
    @callback(
        [Output('filter-a-status', 'children'),
        Output('filter-b-status', 'children'),
        Output('filter-c-status', 'children'),
        Output('filter-d-status', 'children')],
        [Input('integrated-selected-traits-store', 'data'),
        Input('integrated-secondary-filter-enabled', 'data'),
        Input('integrated-trait-groups-store', 'data'),
        Input('integrated-pvalue-cutoff', 'value'),
        Input('integrated-maf-enabled', 'value'),
        Input('integrated-maf-cutoff', 'value'),
        Input('integrated-unique-position-enabled', 'value')],
        prevent_initial_call=True
    )
    def update_filter_tab_status(selected_traits, secondary_enabled, trait_groups, 
                                pvalue_cutoff, maf_enabled, maf_cutoff, unique_enabled):
        """필터 탭의 상태 표시 업데이트"""
        
        # Filter A status (GWAS Related Filter) - Fixed logic to handle grouped traits
        trait_count = len(selected_traits) if selected_traits else 0
        group_count = len(trait_groups.get('groups', {})) if trait_groups else 0
        
        # Check if any groups have traits
        has_grouped_traits = False
        if trait_groups and trait_groups.get('groups'):
            for group_data in trait_groups['groups'].values():
                if group_data is not None and len(group_data) > 0:
                    has_grouped_traits = True
                    break
        
        # Determine status based on traits and groups
        if secondary_enabled and has_grouped_traits:
            a_status = "Trait Grouping"
        elif trait_count > 0:
            a_status = "Trait Selection"
        elif has_grouped_traits:
            # Traits exist in groups but secondary filtering might not be enabled
            a_status = "Trait Grouping"
        else:
            a_status = "Not Used"
        
        # Filter B status (P-value Filter) - 항상 활성화
        b_status = f"P-value ≥ {pvalue_cutoff}" if pvalue_cutoff is not None else "P-value ≥ 5"
        
        # Filter C status (MAF Filter)
        if maf_enabled:
            c_status = f"MAF ≥ {maf_cutoff}"
        else:
            c_status = "Not Used"
        
        # Filter D status (Unique Sample Filter)
        if unique_enabled:
            d_status = "Unique Only"
        else:
            d_status = "Not Used"
        
        return a_status, b_status, c_status, d_status
    '''
    # Filter A detailed info display when collapsed
    
    # Filter Summary Display 콜백 - Updated to show only GWAS Related Filter info
    '''
    @callback(
        Output('filter-summary2', 'children', allow_duplicate=True),
        [
            Input('filter-options-collapse', 'is_open'),  # 필터 collapse 상태
            Input('integrated-subtrait-dropdown', 'value'),
            Input('pvalue-cutoff', 'value'),
            Input('maf-enabled', 'value'),
            Input('maf-cutoff', 'value'),
            Input('snp-occurrence-enabled', 'value'),
            Input('snp-occurrence-min-count', 'value'),
        ],
        prevent_initial_call=True
    )
    def update_filter_summary(is_open,
                            selected_subtraits,
                            pvalue_cutoff, maf_enabled, maf_cutoff,
                            snp_enabled, snp_threshold):
        from dash import html

        # collapse 열려있으면 summary 숨김
        if is_open:
            return ""

        summary_lines = []

        # --- Subtraits ---
        selected_subtraits = selected_subtraits or []
        if selected_subtraits:
            clean_subtraits = [remove_parentheses_numbers(st) for st in selected_subtraits]
            summary_lines.append(
                html.Pre("Subtraits: " + ", ".join(clean_subtraits),
                        style={"margin": 0, "whiteSpace": "pre-wrap"})
            )

        # --- P-value ---
        cutoff_val = float(pvalue_cutoff) if pvalue_cutoff is not None else 5.0
        if cutoff_val != 5.0:
            summary_lines.append(
                html.Pre(f"-log₁₀(P) ≥ {cutoff_val}",
                        style={"margin": 0, "whiteSpace": "pre-wrap"})
            )

        # --- MAF ---
        if maf_enabled:
            summary_lines.append(
                html.Pre(f"MAF ≥ {maf_cutoff}",
                        style={"margin": 0, "whiteSpace": "pre-wrap"})
            )

        # --- SNP Presence ---
        if snp_enabled:
            summary_lines.append(
                html.Pre(f"SNPs in ≥ {snp_threshold} samples",
                        style={"margin": 0, "whiteSpace": "pre-wrap"})
            )

        if not summary_lines:
            return ""

        return html.Div(summary_lines, style={"fontSize": "12px"})
    '''


    
    # Filter Summary Toggle Callback - Improved interaction

     
     # Integrated Subtrait dropdown 옵션 업데이트
    
def create_enhanced_clickable_barplot(summary_df, vcf_samples, selected_traits=None, mode="default", sample_label_dict=None):
    """
    Enhanced barplot generator (Subtrait 구역 구분 포함)
    - mode='default': each / merge (기존)
    - mode='unique': unique mode (샘플별 stacked bar)
    - Subtrait 라벨을 x축에 붙여 표시 (45도, _ → 줄바꿈)
    """
    from dash import html
    import plotly.graph_objects as go
    import plotly.express as px
    import numpy as np
    sample_label_dict = sample_label_dict or {}
    print(sample_label_dict)
    if summary_df is None or summary_df.empty:
        return html.Div("No summary data available", className="text-muted text-center p-3")

    selected_traits = selected_traits or []
    palette = px.colors.qualitative.Set2
    palette2 = px.colors.qualitative.Set3         # Sample용 ✅ 추가

    sample_colors = {s: palette2[i % len(palette2)] for i, s in enumerate(vcf_samples)}

    # 기본 컬럼 확보
    if "Subtrait" not in summary_df.columns:
        summary_df["Subtrait"] = "Unknown"

    # ─────────────────────────────
    # UNIQUE 모드 (stacked bar)
    # ─────────────────────────────
    if mode == "unique":
        fig = go.Figure()
        for s in vcf_samples:
            col = f"{s}_Variant Count" if f"{s}_Variant Count" in summary_df.columns else "Variant Count"
            subdf = summary_df[summary_df.get("Sample", s) == s] if "Sample" in summary_df.columns else summary_df
            if subdf.empty:
                continue
            labels = sample_label_dict.get(s, s) 
            fig.add_bar(
                x=subdf["Trait"],
                y=subdf[col],
                name=labels,
                marker=dict(color=sample_colors.get(s, "#999")),
                text=subdf[col],
                textposition="auto",
                hovertemplate=f"<b>{s}</b><br>%{{x}}: %{{y}}<extra></extra>"
            )
        fig.update_layout(barmode="stack")

    # ─────────────────────────────
    # 기본 (each / merge)
    # ─────────────────────────────
    else:
        is_single = 'Variant Count' in summary_df.columns
        ycol = 'Variant Count' if is_single else 'Total Unique Variants'

        # 색상 설정
        # 색상 설정 (✅ selected_traits 반영)
        colors, opacities = [], []
        for _, row in summary_df.iterrows():
            trait = row.get('Trait', '')
            subtrait = row.get('Subtrait', '')
            try:
                base_color = get_subtrait_color(subtrait)
            except Exception:
                base_color = '#3498db'

            # ✅ 선택된 trait만 강조 (opacity=1)
            if (not selected_traits) or (trait in selected_traits):
                colors.append(base_color)
                opacities.append(1.0)
            else:
                colors.append(base_color)
                opacities.append(0.3)

        # Subtrait별 구간 인덱스 생성
        subtraits = summary_df["Subtrait"].unique()
        x_positions, subtrait_centers = [], []
        group_start = 0
        for st in subtraits:
            subdf = summary_df[summary_df["Subtrait"] == st]
            traits = subdf["Trait"].tolist()
            positions = np.arange(group_start, group_start + len(traits))
            x_positions.extend(positions)
            subtrait_centers.append((positions[0] + positions[-1]) / 2)
            group_start += len(traits) + 2  # 간격 여백

        summary_df = summary_df.copy()
        summary_df["x_pos"] = x_positions[:len(summary_df)]

        fig = go.Figure()

        # Summary bar
        bar_name = "Summary (bar)"
        if is_single and vcf_samples:
            single_sample = vcf_samples[0]
            bar_name = sample_label_dict.get(single_sample, single_sample)  # ✅ 샘플명으로 대체

        fig.add_bar(
            x=summary_df["x_pos"],
            y=summary_df[ycol],
            name=bar_name,
            marker=dict(color=colors, opacity=opacities, line=dict(width=0.5, color="#2c3e50")),
            text=summary_df[ycol],
            textposition="auto",
            textfont=dict(size=11, color="black"),
            customdata=summary_df["Trait"],
            hovertemplate="<b>%{customdata}</b><br>Count: %{y}<extra></extra>",
            showlegend=False  # ✅ 실제 bar는 legend 제외
        )
        summary_color = "#323232"  # 회색 tone for Summary (bar)
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode="markers",
            marker=dict(size=10, color=summary_color, symbol="square"),
            name=bar_name,
            showlegend=True
        ))

        # Sample traces
        if not is_single:
            for s in vcf_samples:
                col = f"{s}_Variant Count"
                if col in summary_df.columns:
                    label = sample_label_dict.get(s, s)  # ✅ dict 적용
                    fig.add_trace(go.Scatter(
                        x=summary_df["x_pos"],
                        y=summary_df[col],
                        mode="lines",
                        name=label,
                        line=dict(width=1),
                        hovertemplate=f"<b>{label}</b><br>%{{customdata}}: %{{y}}<extra></extra>",
                        customdata=summary_df["Trait"]
                    ))

        # Subtrait 배경 + tick 정의
        offset = 0
        subtrait_tickvals = []
        subtrait_ticktext = []
        for st in subtraits:
            subdf = summary_df[summary_df["Subtrait"] == st]
            if subdf.empty:
                continue
            group_end = offset + len(subdf["Trait"].unique()) - 0.5
            center = (offset + group_end) / 2
            subtrait_tickvals.append(center)
            subtrait_ticktext.append(st.replace("_", "<br>"))  # _ → 줄바꿈
            fig.add_vrect(
                x0=offset - 0.5,
                x1=group_end + 1,
                fillcolor="lightgray",
                opacity=0.08,
                layer="below",
                line_width=0,
            )
            offset += len(subdf["Trait"].unique()) + 2

        # ✅ x축 눈금(Subtrait)으로 표시
        fig.update_xaxes(
            tickvals=subtrait_tickvals,
            ticktext=subtrait_ticktext,
            tickangle=45,
            tickfont=dict(size=11, color="black"),
            showgrid=False,
            ticks="outside",
            ticklen=8
        )

    # ─────────────────────────────
    # 공통 레이아웃
    # ─────────────────────────────
    fig.update_layout(
        #title="Trait Barplot (grouped by Subtrait)",
        xaxis_title="Subtrait Name",
        yaxis_title="Variant Count",
        height=400,
        margin=dict(l=60, r=30, t=80, b=120),
        plot_bgcolor="white",
        paper_bgcolor="white",
        hovermode="closest",
        clickmode="event",
        showlegend=True
    )
    fig.update_yaxes(showgrid=True, gridcolor="#e6e6e6", rangemode="tozero")

    return fig





@callback(
    [
        Output('integrated-trait-barplot', 'figure'),
        Output('filter-meta-store', 'data'),
        Output('df-vp-store', 'data'),
    ],
    [
        Input('gwas-selected-samples-store', 'data'),
        Input('integrated-subtrait-dropdown', 'value'),
        Input('integrated-selected-traits-store', 'data'),
        Input('gwas-combo-store', 'data'),
        Input('snp-occurrence-store', 'data'),
        Input('sample-mode-store', 'data'),      # ✅ 현재 모드
        Input('maf-enabled', 'value'),           # ✅ 스위치 입력 추가
    ],
    [
        State('maf-filter-store', 'data'),       # ✅ 필터 상태는 State로 전환
        State('filter-meta-store', 'data'),
        State('gwas-selected-label-store','data')
    ],
    prevent_initial_call=True
)
def update_integrated_trait_bar_container(
    selected_samples, selected_subtraits, selected_traits,
    combo_json, presence_state, mode, maf_enabled,
    maf_filter_state, filter_meta_state,gwas_label
):
    print("🔄 Updating trait bar container...")

    # 기본 방어
    selected_samples = selected_samples or []
    selected_subtraits = selected_subtraits or []
    selected_traits = selected_traits or []
    filter_meta_state = filter_meta_state or {}

    # --- 1️⃣ Validation
    if not selected_samples:
        return html.Div("Please select samples", style={'textAlign': 'center', 'padding': '20px'}), filter_meta_state, None
    if not GWAS_AVAILABLE:
        return html.Div("GWAS module not available", style={'textAlign': 'center', 'padding': '20px'}), filter_meta_state, None
    if not combo_json or combo_json.get('status') != 'ready':
        return html.Div("GWAS data not ready", style={'textAlign': 'center', 'padding': '20px'}), filter_meta_state, None
    
    filter_meta = filter_meta_state.copy()
    # --- 2️⃣ MAF / Presence 필터 설정 (정확한 선언 방식)
    maf_cutoff_default = 0.05
    maf_cutoff = maf_cutoff_default
    maf_active = False

    if maf_enabled:  # ✅ 스위치가 ON인 경우만 활성화
        if maf_filter_state:
            maf_cutoff = maf_filter_state.get("cutoff", maf_cutoff_default)
            maf_active = bool(maf_filter_state.get("enabled", True))
        else:
            maf_cutoff = maf_cutoff_default
            maf_active = True
    else:
        maf_active = False
        maf_cutoff = maf_cutoff_default

    presence_enabled = presence_state.get("enabled", False) if presence_state else False
    presence_threshold = presence_state.get("threshold", None) if presence_state else None

    try:
        # --- 3) 데이터 생성
        summary_df, df_vp2, filter_meta_new = create_sample_summary_data_from_combo(
            combo_json, selected_samples,
            maf_cutoff=maf_cutoff, maf_enabled=maf_enabled,
            presence_threshold=presence_threshold, presence_enabled=presence_enabled
        )

        if summary_df is None or summary_df.empty:
            return html.Div("No summary data", style={'textAlign':'center','padding':'20px'}), filter_meta, None

        # --- 4) Subtrait 필터
        df_plot = summary_df
        df_vp = df_vp2
        if selected_subtraits and 'Subtrait' in summary_df.columns:
            df_plot = summary_df[summary_df['Subtrait'].isin(selected_subtraits)].copy()
            df_vp = df_vp2[df_vp2['Subtrait'].isin(selected_subtraits)].copy()
            if df_plot.empty:
                return html.Div(
                    f"No traits for subtraits: {', '.join(selected_subtraits)}",
                    style={'textAlign':'center','padding':'20px'}
                ), filter_meta, None

                # --- 7) Presence set 누적 저장

            
        # --- 5) 그래프 생성 (모드에 따라 다르게)
        if mode == "unique":
            print("🧩 Rendering UNIQUE barplot (stacked mode)")
            df_plot = df_plot.drop(columns=["Subtrait"], errors="ignore")
            df_plot = df_plot.rename(columns={"Clean_Subtrait": "Subtrait"})
            

            # unique 모드에서는 create_enhanced_clickable_barplot()의 stacked logic을 활용
            graph = create_enhanced_clickable_barplot(
                df_plot,
                selected_samples,
                selected_traits,
                mode="unique"   # ✅ 모드 인자로 전달
                , sample_label_dict=gwas_label
            )

        else:
            print("📊 Rendering regular barplot")
            df_plot = df_plot.drop(columns=["Subtrait"], errors="ignore")
            df_plot = df_plot.rename(columns={"Clean_Subtrait": "Subtrait"})
            graph = create_enhanced_clickable_barplot(
                df_plot,
                selected_samples,
                selected_traits,
                mode="default"
                , sample_label_dict=gwas_label
            )
        
            if presence_enabled and filter_meta_new and "presence" in filter_meta_new:
                presence_info = filter_meta_new["presence"]
                presence_sets = filter_meta.get("presence_sets", [])
                # 중복 방지: threshold 기준으로 unique 처리
                if not any(p["threshold"] == presence_info["threshold"] for p in presence_sets):
                    presence_sets.append(presence_info)
                filter_meta["presence_sets"] = presence_sets

        return graph, filter_meta, df_vp.to_dict('records')

    except Exception as e:
        print(f"❌ Error in bar container callback: {e}")
        return html.Div(f"Error: {e}", style={'color':'red','textAlign':'center'}), filter_meta, None


@callback(
    Output("processed-df-store", "data"),   # ✅ 유일한 출력
    Input("df-vp-store", "data"),           # ✅ GWAS 원본 데이터
    State("sample-mode-store", "data"),     # ✅ 현재 분석 모드 (each / merge / unique)
    prevent_initial_call=True
)
def preprocess_df_vp(df_vp_data, mode):
    import pandas as pd

    # ─────────────────────────────
    # 유효성 검사
    # ─────────────────────────────
    if not df_vp_data:
        print("⚠️ No df_vp_data provided.")
        return {}

    df_vp = pd.DataFrame(df_vp_data)
    base_cols = [
        'Variant ID', 'Chromosome', 'Position', 'Minor Allele', 'MAF',
        'P-value', 'PMID', 'Trait', 'Subtrait', 'Clean_Subtrait'
    ]
    gt_cols = [c for c in df_vp.columns if c.endswith('_GT')]
    samples = [c.replace('_GT', '') for c in gt_cols]

    if not gt_cols:
        print("⚠️ No genotype columns (_GT) found.")
        return {}

    results = {}

    # ─────────────────────────────
    # EACH (각 샘플별 분리)
    # ─────────────────────────────
    each_dict = {}
    for col in gt_cols:
        sample = col.replace('_GT', '')
        mask = df_vp[col].notna()
        sample_df = df_vp.loc[mask, base_cols + [col]].copy()
        each_dict[sample] = sample_df.to_dict('records')
    results["each"] = each_dict

    # ─────────────────────────────
    # MERGE (모든 샘플 통합)
    # ─────────────────────────────
    results["merge"] = df_vp[base_cols + gt_cols].to_dict('records')

    # ─────────────────────────────
    # UNIQUE (샘플 고유 변이)
    # ─────────────────────────────
    if len(gt_cols) > 5:
        results["unique"] = "❌ UNIQUE mode supports ≤5 samples only"
    else:
        uniq_dict = {}
        for col in gt_cols:
            sample = col.replace('_GT', '')
            others = [c for c in gt_cols if c != col]
            mask = (df_vp[col].notna()) & (df_vp[others].isna().all(axis=1))
            uniq_df = df_vp.loc[mask, base_cols + [col]].copy()
            uniq_dict[sample] = uniq_df.to_dict('records')
        results["unique"] = uniq_dict

    # ─────────────────────────────
    # 모드 기반 반환
    # ─────────────────────────────
    if not mode:
        print("ℹ️ No mode set yet — returning all results (each, merge, unique)")
        return results

    # 필요한 모드만 반환
    if mode in results:
        print(f"✅ Returning processed data for mode = '{mode}'")
        return {mode: results[mode]}
    else:
        print(f"⚠️ Mode '{mode}' not found, returning all results.")
        return results

@callback(
    Output("selected-traits-section", "style"),
    Output("selected-traits-summary", "children"),
    Input("copy2selected-traits-store", "data"),
    prevent_initial_call=True
)
def update_selected_traits_summary(selected_traits):
    from dash import html
    selected_traits = selected_traits or []

    if not selected_traits:
        return {"display": "none"}, html.P(
            "No traits selected.",
            style={'fontSize': '12px', 'color': '#777'}
        )

    # ✅ 3개씩 줄바꿈 처리
    trait_chunks = [", ".join(selected_traits[i:i+3]) for i in range(0, len(selected_traits), 3)]
    summary_divs = [
        html.Div(chunk, style={
            "fontFamily": "monospace",
            "color": "#333",
            "marginBottom": "2px",
            "fontSize": "13px",
            "textAlign": "left"
        }) for chunk in trait_chunks
    ]

    box_style = {
        "display": "block",
        "border": "1px solid #ccc",
        "borderRadius": "4px",
        "padding": "8px",
        "width": "100%",
        "margin": "0 auto",
        "backgroundColor": "#f8f9fa",
        "boxShadow": "0 1px 2px rgba(0,0,0,0.05)",
    }

    return box_style, summary_divs

@callback(
    Output("integrated-selected-traits-table", "style"),
    Output("selected-traits-status", "children"),
    Input("selected-traits-toggle", "n_clicks"),
    State("integrated-selected-traits-table", "style"),
    prevent_initial_call=True
)
def toggle_table_visibility(n_clicks, current_style):
    if not current_style:
        current_style = {
            'marginTop': '10px',
            'maxHeight': '0px',
            'overflow': 'hidden',
            'transition': 'max-height 0.3s ease-out'
        }

    if current_style.get("maxHeight") == "0px":
        current_style["maxHeight"] = "600px"
        return current_style, "(opened)"
    else:
        current_style["maxHeight"] = "0px"
        return current_style, "(closed)"

# --- (2) 테이블 내용 업데이트 ---
@callback(
    Output("integrated-selected-traits-table", "children"),
    Input("copy2selected-traits-store", "data"),
    prevent_initial_call=True
)
def update_selected_traits_table(selected_traits):
    from dash import html
    selected_traits = selected_traits or []

    if not selected_traits:
        return html.P("No traits selected.", style={'fontSize': '12px', 'color': '#777'})

    try:
        return create_selected_traits_table(selected_traits, TRAIT_INFO_DF)
    except Exception as e:
        return html.P(f"Error creating table: {str(e)}", style={'color': 'red'})

    '''
    # 통합 테이블 업데이트 (그룹 컬럼 포함)
    @callback(
        [Output('gwas-locus-table', 'columns', allow_duplicate=True),
         Output('gwas-locus-table', 'data', allow_duplicate=True)],
        [Input('integrated-selected-traits-store', 'data'),
         Input('integrated-secondary-filter-enabled', 'data'),
         Input('integrated-trait-groups-store', 'data'),
         Input('copy2gwas-selected-samples-store', 'data')],
        prevent_initial_call=True
    )
    def update_integrated_locus_table(selected_traits, secondary_enabled, trait_groups, selected_samples):
        """통합 필터 시스템에서 테이블 업데이트 (그룹 컬럼 포함)"""
        
        # Check if traits are available (either selected or grouped)
        has_traits = False
        all_traits = []
        
        if selected_traits and len(selected_traits) > 0:
            has_traits = True
            all_traits = selected_traits.copy()
        
        # Add grouped traits to the mix
        if trait_groups and trait_groups.get('groups'):
            for group_data in trait_groups['groups'].values():
                if group_data is not None and len(group_data) > 0:
                    has_traits = True
                    all_traits.extend(group_data)
        
        # Remove duplicates
        all_traits = list(set(all_traits))
        
        if not selected_samples or not has_traits or not GWAS_AVAILABLE:
            print(f"❌ Table update skipped: samples={bool(selected_samples)}, has_traits={has_traits}, GWAS_AVAILABLE={GWAS_AVAILABLE}")
            return [], []
        
        try:
            gwas_result = get_gwas_data_for_samples_optimized(selected_samples)
            if not gwas_result:
                return [], []
            
            gwas_module = gwas_result['gwas_module']
            
            # 기본 테이블 생성 - use all available traits
            base_columns, base_data = gwas_module.update_locus_table_data(
                all_traits, selected_samples, None, {}, []
            )
            
            print(f"✅ Table data generated: {len(base_data)} rows for {len(all_traits)} traits")
            
            # 2차 필터링이 활성화된 경우 그룹 컬럼 추가
            if secondary_enabled and trait_groups and trait_groups.get('groups'):
                # 그룹 정보 파싱
                groups = trait_groups['groups']
                group_trait_map = {}  # trait -> group_number 매핑
                
                for group_name, traits_in_group in groups.items():
                    if traits_in_group is not None and len(traits_in_group) > 0:
                        group_number = group_name.replace('Group ', '').replace('group', '').strip()
                        for trait in traits_in_group:
                            group_trait_map[trait] = f"group{group_number}"
                
                # 그룹 컬럼 추가
                if group_trait_map:
                    # 새로운 컬럼 정의
                    enhanced_columns = base_columns.copy()
                    
                    # 사용된 그룹들 확인
                    used_groups = set(group_trait_map.values())
                    for group_col in sorted(used_groups):
                        enhanced_columns.append({
                            'name': group_col.capitalize(),
                            'id': group_col,
                            'type': 'text'
                        })
                    
                    # 데이터에 그룹 정보 추가
                    enhanced_data = []
                    for row in base_data:
                        enhanced_row = row.copy()
                        trait = row.get('Trait', '')
                        
                        # 각 그룹에 대해 체크
                        for group_col in sorted(used_groups):
                            if group_trait_map.get(trait) == group_col:
                                enhanced_row[group_col] = '✓'
                            else:
                                enhanced_row[group_col] = ''
                        
                        enhanced_data.append(enhanced_row)
                    
                    print(f"   ✅ 통합 테이블 생성 (그룹 포함): {len(enhanced_columns)}열, {len(enhanced_data)}행")
                    return enhanced_columns, enhanced_data
            
            # 2차 필터링이 비활성화된 경우 기본 테이블 반환
            print(f"   ✅ 통합 테이블 생성 (기본): {len(base_columns)}열, {len(base_data)}행")
            return base_columns, base_data
            
        except Exception as e:
            print(f"❌ 통합 테이블 업데이트 오류: {e}")
            import traceback
            traceback.print_exc()
            return [], []
    '''

    '''
    # 그룹에서 trait 제거 (모든 그룹 지원) - Selected Traits에 다시 추가
    @callback(
        [Output('integrated-trait-groups-store', 'data', allow_duplicate=True),
         Output('integrated-trait-group-1', 'children', allow_duplicate=True),
         Output('integrated-trait-group-2', 'children', allow_duplicate=True),
         Output('integrated-dynamic-groups-container', 'children', allow_duplicate=True),
         Output('integrated-selected-traits-store', 'data', allow_duplicate=True)],
        Input({'type': 'integrated-group-chip', 'trait': ALL, 'group': ALL}, 'n_clicks'),
        [State('integrated-trait-groups-store', 'data'),
         State('integrated-selected-traits-store', 'data'),
         State('integrated-dynamic-groups-container', 'children')],
        prevent_initial_call=True
    )
    def remove_trait_from_group(chip_clicks, current_groups, selected_traits, current_dynamic_ui):
        """그룹에서 trait 제거하고 Selected Traits에 다시 추가"""
        ctx = dash.callback_context
        
        if not ctx.triggered or not any(chip_clicks):
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update
        
        # 클릭된 chip 확인
        triggered_id = ctx.triggered[0]['prop_id']
        if 'integrated-group-chip' not in triggered_id:
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update
        
        import json
        chip_info = json.loads(triggered_id.split('.')[0])
        trait_to_remove = chip_info['trait']
        
        print(f"🗑️ Removing trait '{trait_to_remove}' from groups")
        
        # Selected Traits에 다시 추가
        updated_selected_traits = selected_traits.copy() if selected_traits else []
        if trait_to_remove not in updated_selected_traits:
            updated_selected_traits.append(trait_to_remove)
            print(f"   Trait '{trait_to_remove}' restored to Selected Traits")
        
        # 그룹 데이터에서 trait 제거
        updated_groups = current_groups.copy() if current_groups else {'groups': {}, 'max_groups': 9}
        if 'groups' not in updated_groups:
            updated_groups['groups'] = {}
        
        for group_name, traits in updated_groups['groups'].items():
            if traits is not None and len(traits) > 0 and trait_to_remove in traits:
                traits.remove(trait_to_remove)
                print(f"   ✅ Removed '{trait_to_remove}' from {group_name}")
                break
        
        # 그룹 UI 업데이트
        def create_group_chips(traits, color):
            if not traits:
                return [html.P("Click 'add trait+' to add traits here...", style={'color': '#999', 'text-align': 'center', 'margin': '0', 'font-style': 'italic'})]
            
            chips = []
            for trait in traits:
                chip = dbc.Button([
                    trait + " ",
                    html.Span("×", style={'margin-left': '5px', 'font-weight': 'bold'})
                ], 
                id={'type': 'integrated-group-chip', 'trait': trait, 'group': color},
                size='sm', 
                color='light',
                style={
                    'margin': '2px',
                    'border': f'1px solid {color}',
                    'background-color': 'white'
                })
                chips.append(chip)
            return chips
        
        group1_ui = create_group_chips(updated_groups['groups'].get('Group 1', []), px.colors.qualitative.Pastel1[0])
        group2_ui = create_group_chips(updated_groups['groups'].get('Group 2', []), px.colors.qualitative.Pastel1[1])
        
        # 동적 그룹들 UI 업데이트 (Group 3부터)
        updated_dynamic_ui = []
        
        # 기존 그룹들을 모두 재생성 (Group 3부터)
        if updated_groups and 'groups' in updated_groups:
            groups = updated_groups['groups']
            for group_name, group_traits in groups.items():
                if group_name.startswith('Group ') and group_name not in ['Group 1', 'Group 2']:
                    try:
                        group_num = int(group_name.split(' ')[1])
                        
                        # 그룹 색상 결정
                        group_color_idx = group_num - 1
                        if group_color_idx < len(px.colors.qualitative.Pastel1):
                            group_color = px.colors.qualitative.Pastel1[group_color_idx]
                        else:
                            group_color = '#999'
                        
                        # Use consistent colored circle symbols
                        
                        # 그룹 UI 생성
                        new_group_div = html.Div([
                            html.Div([
                                 dbc.Button([
                                    html.Span("●", style={'color': group_color, 'font-size': '14px', 'margin-right': '8px'}),
                                    group_name
                                 ], 
                                           id={'type': 'integrated-add-to-group', 'group': group_num}, 
                                           size="sm", 
                                           color="info",
                                           outline=True,
                                           style={'font-weight': 'bold', 'width': '100%', 'margin-bottom': '5px'}),
                                 html.P("Click to add traits →", style={'font-size': '0.8em', 'color': '#666', 'margin': '0 0 10px 0', 'text-align': 'center'})
                             ]),
                            html.Div(
                                id=f'integrated-trait-group-{group_num}',
                                style={
                                    'min-height': '60px',
                                    'border': f'2px dashed {group_color}',
                                    'border-radius': '8px',
                                    'padding': '15px',
                                    'margin-bottom': '15px',
                                    'background-color': '#f8f9fa'
                                },
                                children=create_group_chips(group_traits, group_color)
                            )
                        ], style={'width': '30%', 'display': 'inline-block', 'margin-right': '3%', 'vertical-align': 'top'})
                        
                        updated_dynamic_ui.append(new_group_div)
                    except:
                        continue
        
        return updated_groups, group1_ui, group2_ui, updated_dynamic_ui, updated_selected_traits
    
    # 모달 열기/닫기 (모든 Add 버튼 동적 처리)
    @callback(
        [Output('integrated-trait-selection-modal', 'is_open', allow_duplicate=True),
         Output('integrated-selected-group-info', 'children', allow_duplicate=True),
         Output('integrated-trait-selection-dropdown', 'options', allow_duplicate=True),
         Output('integrated-trait-selection-dropdown', 'value', allow_duplicate=True),
         Output('integrated-group-selection-dropdown', 'options', allow_duplicate=True),
         Output('integrated-group-selection-dropdown', 'value', allow_duplicate=True)],
        [Input({'type': 'integrated-add-to-group', 'group': ALL}, 'n_clicks'),  # 모든 동적 Add 버튼
         Input('integrated-add-to-group-1', 'n_clicks'),
         Input('integrated-add-to-group-2', 'n_clicks'),
         Input('integrated-trait-modal-cancel', 'n_clicks'),
         Input('integrated-trait-modal-confirm', 'n_clicks')],
        [State('integrated-selected-traits-store', 'data'),
         State('integrated-trait-groups-store', 'data')],
        prevent_initial_call=True
    )
    def handle_trait_add_modal(dynamic_clicks, group1_clicks, group2_clicks, cancel_clicks, confirm_clicks,
                              selected_traits, trait_groups):
        """Trait 추가 모달 처리 (모든 Add 버튼 포함 - 동적 그룹 지원)"""
        ctx = dash.callback_context
        
        if not ctx.triggered:
            return False, '', [], None, [], None
        
        # 실제 클릭이 발생했는지 확인 (n_clicks > 0)
        trigger_value = ctx.triggered[0]['value']
        if not trigger_value or trigger_value == 0:
            return False, '', [], None, [], None
        
        triggered_prop = ctx.triggered[0]['prop_id']
        triggered_id = triggered_prop.split('.')[0]
        
        # 동적 Add 버튼 (패턴 매칭) 또는 기본 Add 버튼 처리
        group_number = None
        
        # 동적 그룹 버튼 확인 (패턴 매칭)
        if '"type":"integrated-add-to-group"' in triggered_id:
            try:
                import json
                button_info = json.loads(triggered_id)
                group_number = str(button_info['group'])
            except:
                pass
        
        # 기본 그룹 버튼 확인
        elif 'integrated-add-to-group-' in triggered_id:
            group_number = triggered_id.split('-')[-1]  # '1', '2'
        
        if group_number:
            if not selected_traits:
                error_info = html.P("No traits selected. Please select traits first.", style={'color': 'red'})
                return False, error_info, [], None, [], None
            
            group_name = f"Group {group_number}"
            
            # 선택한 그룹에 이미 있는 trait들만 제외 (다른 그룹에 있는 것은 허용)
            groups = trait_groups.get('groups', {}) if trait_groups else {}
            current_group_traits = groups.get(group_name, [])
            
            available_traits = [trait for trait in selected_traits if trait not in current_group_traits]
            
            if not available_traits:
                warning_info = html.P(f"All selected traits are already in {group_name}.", style={'color': 'orange'})
                return False, warning_info, [], None, [], None
            
            # 모달 열기
            trait_options = [{'label': trait, 'value': trait} for trait in available_traits]
            
            # Use consistent colored circle symbols for groups
            group_num = int(group_number)
            group_color = px.colors.qualitative.Pastel1[group_num - 1] if group_num <= len(px.colors.qualitative.Pastel1) else px.colors.qualitative.Pastel1[0]
            
            # 사용 가능한 그룹 옵션들 동적 생성
            group_options = [
                {'label': 'Group 1', 'value': 'Group 1'},
                {'label': 'Group 2', 'value': 'Group 2'}
            ]
            
            # 기존에 생성된 그룹들도 추가
            if trait_groups and 'groups' in trait_groups:
                for existing_group in trait_groups['groups'].keys():
                    if existing_group not in ['Group 1', 'Group 2']:
                        try:
                            group_num = int(existing_group.split(' ')[1])
                            emoji_existing = emoji_map.get(group_num, '🔘')
                            group_option = {'label': f'{emoji_existing} {existing_group}', 'value': existing_group}
                            if group_option not in group_options:
                                group_options.append(group_option)
                        except:
                            continue
            
            # Create colored circle symbol based on group
            group_num = int(group_name.split(' ')[1]) if 'Group ' in group_name else 1
            group_color = px.colors.qualitative.Pastel1[group_num - 1] if group_num <= len(px.colors.qualitative.Pastel1) else px.colors.qualitative.Pastel1[0]
            
            group_info = html.Div([
                html.P([
                    "Add trait to ",
                    html.Span("●", style={'color': group_color, 'font-size': '14px', 'margin-right': '5px'}),
                    f"{group_name}:"
                ], style={'font-weight': 'bold', 'color': 'blue', 'margin': '0'}),
                html.P("Select a trait and confirm.", 
                       style={'color': '#666', 'font-size': '0.9em', 'margin': '5px 0 0 0'})
            ])
            
            return True, group_info, trait_options, None, group_options, group_name
        
        # 취소 버튼
        elif triggered_id == 'integrated-trait-modal-cancel':
            return False, '', [], None, [], None
        
        # 확인 버튼 - 실제 trait 추가는 다른 callback에서 처리
        elif triggered_id == 'integrated-trait-modal-confirm':
            return False, '', [], None, [], None
        
        return False, '', [], None, [], None
    
    # Trait을 그룹에 추가 (모든 그룹 지원) - Selected Traits에서 제거
    @callback(
        [Output('integrated-trait-groups-store', 'data', allow_duplicate=True),
         Output('integrated-trait-group-1', 'children', allow_duplicate=True),
         Output('integrated-trait-group-2', 'children', allow_duplicate=True),
         Output('integrated-dynamic-groups-container', 'children', allow_duplicate=True),
         Output('integrated-selected-traits-store', 'data', allow_duplicate=True)],
        Input('integrated-trait-modal-confirm', 'n_clicks'),
        [State('integrated-trait-selection-dropdown', 'value'),
         State('integrated-group-selection-dropdown', 'value'),
         State('integrated-trait-groups-store', 'data'),
         State('integrated-dynamic-groups-container', 'children'),
         State('integrated-selected-traits-store', 'data')],
        prevent_initial_call=True
    )
    def add_trait_to_group(confirm_clicks, selected_trait, target_group, current_groups, current_dynamic_ui, selected_traits):
        """선택된 trait을 target 그룹에 추가하고 Selected Traits에서 제거"""
        if not selected_trait or not target_group:
            # 변화 없음
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update
        
        # 그룹 데이터 업데이트
        updated_groups = current_groups.copy() if current_groups else {'groups': {}, 'max_groups': 9}
        if 'groups' not in updated_groups:
            updated_groups['groups'] = {}
        
        # 선택된 trait을 target 그룹에 추가
        if target_group not in updated_groups['groups']:
            updated_groups['groups'][target_group] = []
        
        if selected_trait not in updated_groups['groups'][target_group]:
            updated_groups['groups'][target_group].append(selected_trait)
        
        print(f"✅ Trait '{selected_trait}' added to {target_group}")
        
        # Selected Traits에서 해당 trait 제거
        updated_selected_traits = selected_traits.copy() if selected_traits else []
        if selected_trait in updated_selected_traits:
            updated_selected_traits.remove(selected_trait)
            print(f"   Trait '{selected_trait}' removed from Selected Traits")
        
        # 그룹 UI 업데이트 함수
        def create_group_chips(traits, color):
            if not traits:
                return [html.P("Click 'add trait+' to add traits here...", style={'color': '#999', 'text-align': 'center', 'margin': '0', 'font-style': 'italic'})]
            
            chips = []
            for trait in traits:
                chip = dbc.Button([
                    trait + " ",
                    html.Span("×", style={'margin-left': '5px', 'font-weight': 'bold'})
                ], 
                id={'type': 'integrated-group-chip', 'trait': trait, 'group': color},
                size='sm', 
                color='light',
                style={
                    'margin': '2px',
                    'border': f'1px solid {color}',
                    'background-color': 'white'
                })
                chips.append(chip)
            return chips
        
        # Group 1과 Group 2 UI 업데이트
        group1_ui = create_group_chips(updated_groups['groups'].get('Group 1', []), px.colors.qualitative.Pastel1[0])
        group2_ui = create_group_chips(updated_groups['groups'].get('Group 2', []), px.colors.qualitative.Pastel1[1])
        
        # 동적 그룹들 UI 업데이트 (더 안정적인 방법)
        updated_dynamic_ui = []
        
        # 기존 그룹들을 모두 재생성 (Group 3부터)
        if updated_groups and 'groups' in updated_groups:
            groups = updated_groups['groups']
            for group_name, group_traits in groups.items():
                if group_name.startswith('Group ') and group_name not in ['Group 1', 'Group 2']:
                    try:
                        group_num = int(group_name.split(' ')[1])
                        
                        # 그룹 색상 결정
                        group_color_idx = group_num - 1
                        if group_color_idx < len(px.colors.qualitative.Pastel1):
                            group_color = px.colors.qualitative.Pastel1[group_color_idx]
                        else:
                            group_color = '#999'
                        
                        # Use consistent colored circle symbols
                        
                        # 그룹 UI 생성
                        new_group_div = html.Div([
                            html.Div([
                                 dbc.Button([
                                    html.Span("●", style={'color': group_color, 'font-size': '14px', 'margin-right': '8px'}),
                                    group_name
                                 ], 
                                           id={'type': 'integrated-add-to-group', 'group': group_num}, 
                                           size="sm", 
                                           color="info",
                                           outline=True,
                                           style={'font-weight': 'bold', 'width': '100%', 'margin-bottom': '5px'}),
                                 html.P("Click to add traits →", style={'font-size': '0.8em', 'color': '#666', 'margin': '0 0 10px 0', 'text-align': 'center'})
                             ]),
                            html.Div(
                                id=f'integrated-trait-group-{group_num}',
                                style={
                                    'min-height': '60px',
                                    'border': f'2px dashed {group_color}',
                                    'border-radius': '8px',
                                    'padding': '15px',
                                    'margin-bottom': '15px',
                                    'background-color': '#f8f9fa'
                                },
                                children=create_group_chips(group_traits, group_color)
                            )
                        ], style={'width': '30%', 'display': 'inline-block', 'margin-right': '3%', 'vertical-align': 'top'})
                        
                        updated_dynamic_ui.append(new_group_div)
                    except:
                        continue
        
        return updated_groups, group1_ui, group2_ui, updated_dynamic_ui, updated_selected_traits
    
    # Trait Pool 업데이트 (선택된 trait을 개별 버튼으로 표시)
    @callback(
        Output('integrated-trait-pool', 'children'),
        [Input('integrated-selected-traits-store', 'data'),
         Input('integrated-trait-groups-store', 'data')]
    )
    def update_trait_pool(selected_traits, trait_groups):
        """선택된 trait들을 pool에 버튼으로 표시 (그룹에 없는 trait만)"""
        if not selected_traits:
            return [html.P("No traits selected for grouping", 
                          style={'color': '#666', 'margin': '0', 'font-style': 'italic'})]
        
        # 이미 그룹에 할당된 trait들 확인
        groups = trait_groups.get('groups', {}) if trait_groups else {}
        assigned_traits = set()
        for group_traits in groups.values():
            if group_traits:
                assigned_traits.update(group_traits)
        
        # 할당되지 않은 trait들만 pool에 표시
        unassigned_traits = [trait for trait in selected_traits if trait not in assigned_traits]
        
        if not unassigned_traits:
            return [html.P("All traits are assigned to groups", 
                          style={'color': '#28a745', 'margin': '0', 'font-style': 'italic'})]
        
        # 각 trait을 개별 버튼으로 생성
        trait_buttons = []
        for trait in unassigned_traits:
            button = dbc.Button(
                trait,
                id={'type': 'integrated-pool-trait', 'trait': trait},
                size='sm',
                color='light',
                outline=True,
                style={
                    'margin': '3px',
                    'border': '1px solid #6c757d',
                    'background-color': '#ffffff'
                }
            )
            trait_buttons.append(button)
        
        return trait_buttons
    
    # Pool에서 trait 클릭 시 그룹 선택 모달 열기
    @callback(
        [Output('integrated-trait-selection-modal', 'is_open', allow_duplicate=True),
         Output('integrated-selected-group-info', 'children', allow_duplicate=True),
         Output('integrated-trait-selection-dropdown', 'options', allow_duplicate=True),
         Output('integrated-trait-selection-dropdown', 'value', allow_duplicate=True),
         Output('integrated-group-selection-dropdown', 'value', allow_duplicate=True)],
        Input({'type': 'integrated-pool-trait', 'trait': ALL}, 'n_clicks'),
        prevent_initial_call=True
    )
    def handle_pool_trait_click(trait_clicks):
        """Pool의 trait 버튼 클릭 시 어느 그룹에 추가할지 선택"""
        ctx = dash.callback_context
        
        if not ctx.triggered or not any(trait_clicks):
            return False, '', [], None, None
        
        # 클릭된 trait 확인
        triggered_id = ctx.triggered[0]['prop_id']
        if 'integrated-pool-trait' not in triggered_id:
            return False, '', [], None, None
        
        import json
        trait_info = json.loads(triggered_id.split('.')[0])
        clicked_trait = trait_info['trait']
        
        # 자동으로 선택된 trait으로 설정
        trait_options = [{'label': clicked_trait, 'value': clicked_trait}]
        
        group_info = html.Div([
            html.P(f"Assign trait '{clicked_trait}' to a group:", 
                   style={'font-weight': 'bold', 'color': 'blue', 'margin': '0'}),
            html.P("Select the target group below and click 'Add to Group'.", 
                   style={'color': '#666', 'font-size': '0.9em', 'margin': '5px 0 0 0'})
        ])
        
        return True, group_info, trait_options, clicked_trait, None
    
    # Add New Group 기능 구현
    @callback(
        [Output('integrated-dynamic-groups-container', 'children'),
         Output('integrated-trait-groups-store', 'data', allow_duplicate=True)],
        Input('integrated-add-new-group-btn', 'n_clicks'),
        [State('integrated-dynamic-groups-container', 'children'),
         State('integrated-trait-groups-store', 'data')],
        prevent_initial_call=True
    )
    def add_new_group(n_clicks, current_groups_ui, current_groups_data):
        """새로운 그룹 추가 (Group 3부터 시작)"""
        if not n_clicks:
            return current_groups_ui, current_groups_data
        
        # 현재 그룹 수 확인
        current_data = current_groups_data if current_groups_data else {'groups': {}, 'max_groups': 9}
        existing_groups = current_data.get('groups', {})
        
        # 다음 그룹 번호 계산 (Group 3부터 시작)
        max_group_num = 2  # 기본 Group 1, 2
        for group_name in existing_groups.keys():
            if group_name.startswith('Group '):
                try:
                    group_num = int(group_name.split(' ')[1])
                    max_group_num = max(max_group_num, group_num)
                except:
                    continue
        
        new_group_num = max_group_num + 1
        new_group_name = f"Group {new_group_num}"
        
        # 최대 9개 그룹 제한
        if new_group_num > 9:
            return current_groups_ui, current_groups_data
        
        # 새 그룹 UI 생성
        group_color = px.colors.qualitative.Pastel1[new_group_num - 1] if new_group_num <= len(px.colors.qualitative.Pastel1) else '#999'
        
        new_group_div = html.Div([
            html.Div([
                 dbc.Button([
                    html.Span("●", style={'color': group_color, 'font-size': '14px', 'margin-right': '8px'}),
                    new_group_name
                 ], 
                           id={'type': 'integrated-add-to-group', 'group': new_group_num}, 
                           size="sm", 
                           color="info",
                           outline=True,
                           style={'font-weight': 'bold', 'width': '100%', 'margin-bottom': '5px'}),
                 html.P("Click to add traits →", style={'font-size': '0.8em', 'color': '#666', 'margin': '0 0 10px 0', 'text-align': 'center'})
             ]),
            html.Div(
                id=f'integrated-trait-group-{new_group_num}',
                style={
                    'min-height': '60px',
                    'border': f'2px dashed {group_color}',
                    'border-radius': '8px',
                    'padding': '15px',
                    'margin-bottom': '15px',
                    'background-color': '#f8f9fa'
                },
                children=[
                    html.P("Click 'add trait+' to add traits here...", style={
                        'color': '#999', 
                        'text-align': 'center',
                        'margin': '0',
                        'font-style': 'italic'
                    })
                ]
            )
        ], style={'width': '30%', 'display': 'inline-block', 'margin-right': '3%', 'vertical-align': 'top'})
        
        # UI 업데이트
        updated_groups_ui = current_groups_ui if current_groups_ui else []
        updated_groups_ui.append(new_group_div)
        
        # 데이터 업데이트 (빈 그룹으로 추가)
        updated_data = current_data.copy()
        updated_data['groups'][new_group_name] = []
        
        print(f"✅ New group added: {new_group_name}")
        
        return updated_groups_ui, updated_data
    '''
    
    # 통합 필터 → Scatter Plot 연동 (2차 필터링 시 그룹 색상 적용)
    '''
    @callback(
        Output('gwas-sample-scatter', 'figure', allow_duplicate=True),
        [Input('integrated-selected-traits-store', 'data'),
         Input('integrated-secondary-filter-enabled', 'data'),
         Input('integrated-trait-groups-store', 'data'),
         Input('integrated-filter-b-state', 'data'),
         Input('integrated-filter-c-state', 'data'),
         Input('integrated-filter-d-state', 'data'),
         Input('copy2gwas-selected-samples-store', 'data')],  # State에서 Input으로 변경
        prevent_initial_call=True
    )
    def update_integrated_scatter(selected_traits, secondary_enabled, trait_groups, 
                                integrated_filter_b, integrated_filter_c, integrated_filter_d, selected_samples):
        """통합 필터에서 선택된 traits를 scatter plot에 반영"""
        print(f"🔍 update_integrated_scatter 호출됨: selected_samples={selected_samples}")

        # Check if traits are available (either selected or grouped)
        has_traits = False
        all_traits = []
        
        if selected_traits and len(selected_traits) > 0:
            has_traits = True
            all_traits = selected_traits.copy()
        
        # Add grouped traits to the mix
        if trait_groups and trait_groups.get('groups'):
            for group_data in trait_groups['groups'].values():
                if group_data is not None and len(group_data) > 0:
                    has_traits = True
                    all_traits.extend(group_data)
        
        # Remove duplicates
        all_traits = list(set(all_traits))
        print(f"🔍 update_integrated_scatter 호출됨: selected_samples={selected_samples} (타입: {type(selected_samples)})")
        if not selected_samples or not GWAS_AVAILABLE:
            print(f"❌ 통합 Scatter update skipped: samples={bool(selected_samples)}, GWAS_AVAILABLE={GWAS_AVAILABLE}")
            return go.Figure().add_annotation(
                text="샘플을 선택해주세요" if not selected_samples else "GWAS 모듈이 사용할 수 없습니다",
                xref="paper", yref="paper",
                x=0.5, y=0.5, xanchor='center', yanchor='middle',
                showarrow=False, font=dict(size=16, color="gray")
            )
        
        try:
            gwas_result = get_gwas_data_for_samples_optimized(selected_samples)
            if not gwas_result:
                return go.Figure().add_annotation(
                    text="GWAS 데이터를 로드할 수 없습니다",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5, xanchor='center', yanchor='middle',
                    showarrow=False, font=dict(size=16, color="red")
                )
            
            # Apply integrated filter states
            filter_states = {}
            pvalue_cutoff = 5  # default
            show_unique_only = False
            
            # Integrated Filter B (P-value cutoff)
            if integrated_filter_b and integrated_filter_b.get('enabled', False):
                pvalue_cutoff = integrated_filter_b.get('pvalue_cutoff', 5)
                print(f"🎯 Integrated P-value cutoff applied: {pvalue_cutoff}")
            
            # Integrated Filter C (MAF cutoff)
            if integrated_filter_c and integrated_filter_c.get('enabled', False):
                filter_states['maf'] = integrated_filter_c.get('maf_cutoff', 0.05)
                filter_states['maf_enabled'] = True  # enabled 상태 전달
                print(f"🎯 Integrated MAF cutoff applied: {filter_states['maf']}")
            elif integrated_filter_c:
                # MAF filter가 정의되어 있지만 비활성화된 경우
                filter_states['maf_enabled'] = False
                print(f"🎯 Integrated MAF filter disabled")
            
            # Integrated Filter D (Unique position)
            if integrated_filter_d and integrated_filter_d.get('enabled', False):
                show_unique_only = True
                print(f"🎯 Integrated unique position filter applied")
            
            # 2차 필터링이 활성화되고 그룹이 있는 경우 그룹 색상 사용
            if secondary_enabled and trait_groups and trait_groups.get('groups'):
                print("🎨 2차 필터링 활성화: 그룹 색상으로 scatter plot 업데이트")
                
                # 그룹 정보 파싱
                groups = trait_groups['groups']
                trait_to_group = {}
                group_colors = {}
                
                for i, (group_name, traits_in_group) in enumerate(groups.items()):
                    if traits_in_group is not None and len(traits_in_group) > 0:
                        group_color = px.colors.qualitative.Pastel1[i % len(px.colors.qualitative.Pastel1)]
                        group_colors[group_name] = group_color
                        for trait in traits_in_group:
                            trait_to_group[trait] = (group_name, group_color)
                
                # Enhanced scatter plot with group colors - use all available traits
                scatter_fig = create_enhanced_gwas_scatter_plot(
                    gwas_result, selected_samples, None, all_traits,
                    filter_states, [], pvalue_cutoff, show_unique_only,
                    use_group_colors=True,
                    trait_to_group=trait_to_group,
                    hide_vals=None,
                    sample_visible={sample: True for sample in selected_samples},
                    subtrait_visible={}
                )
                
                return scatter_fig
            else:
                print("🎨 1차 필터링만 활성화: 기본 subtrait 색상으로 scatter plot 업데이트")
                
                # 기본 scatter plot (subtrait 색상) - use all available traits  
                scatter_fig = create_enhanced_gwas_scatter_plot(
                    gwas_result, selected_samples, None, all_traits,
                    filter_states, [], pvalue_cutoff, show_unique_only,
                    hide_vals=None,
                    sample_visible={sample: True for sample in selected_samples},
                    subtrait_visible={}
                )
                
                return scatter_fig
                
        except Exception as e:
            print(f"❌ 통합 Scatter plot 업데이트 오류: {e}")
            import traceback
            traceback.print_exc()
            return go.Figure().add_annotation(
                text=f"Scatter plot 업데이트 중 오류: {str(e)}",
                xref="paper", yref="paper",
                x=0.5, y=0.5, xanchor='center', yanchor='middle',
                showarrow=False, font=dict(size=16, color="red")
            )
    
    # 통합 필터 B (P-value) 상태 업데이트 - 항상 활성화
    @callback(
        Output('integrated-filter-b-state', 'data'),
        [Input('integrated-pvalue-cutoff', 'value')],
        prevent_initial_call=True
    )
    def update_integrated_filter_b_state(pvalue_cutoff):
        """P-value 필터 값을 store에 반영 - 항상 활성화됨"""
        print(f"🎯 Integrated Filter B 상태 업데이트: 항상 활성화, cutoff={pvalue_cutoff}")
        return {
            'enabled': True,  # 항상 활성화
            'pvalue_cutoff': float(pvalue_cutoff) if pvalue_cutoff is not None else 5.0
        }
    '''
    # 통합 필터 C (MAF) 상태 업데이트
    '''
    @callback(
        Output('integrated-filter-c-state', 'data'),
        [Input('integrated-maf-enabled', 'value'),
         Input('integrated-maf-cutoff', 'value')],
        [State('integrated-filter-c-collapse', 'is_open')],  # collapse 상태 확인
        prevent_initial_call=True
    )
    def update_integrated_filter_c_state(maf_enabled, maf_cutoff, collapse_is_open):
        """MAF 필터 활성화 상태와 cutoff 값을 store에 반영"""
        # 🔥 핵심 수정: Enable MAF filtering이 비활성화되어 있을 때는 threshold 변경을 무시
        ctx = dash.callback_context
        if ctx.triggered:
            trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
            if trigger_id == 'integrated-maf-cutoff' and not maf_enabled:
                print(f"   🔍 MAF cutoff 변경 무시: MAF filtering 비활성화 상태 (threshold={maf_cutoff})")
                # Enable MAF filtering이 비활성화되어 있을 때는 이전 상태 유지
                raise PreventUpdate
        
        print(f"🎯 Integrated Filter C 상태 업데이트: enabled={maf_enabled}, cutoff={maf_cutoff}, collapse_open={collapse_is_open}")
        
        return {
            'enabled': bool(maf_enabled),
            'maf_cutoff': float(maf_cutoff) if maf_cutoff is not None else 0.05
        }
    
    # 통합 필터 D (Unique Position) 상태 업데이트
    @callback(
        Output('integrated-filter-d-state', 'data', allow_duplicate=True),
        [Input('integrated-unique-position-enabled', 'value')],
        prevent_initial_call=True
    )
    def update_integrated_filter_d_state(unique_enabled):
        """Unique position 필터 활성화 상태를 store에 반영"""
        return {
            'enabled': bool(unique_enabled)
        }
    
    # OLD Unique Position 필터 UI 제어 - 새로운 integrated filter로 대체됨

    
    # Filter D 탭 비활성화 (샘플 수 부족시)
    
    @callback(
        [Output('integrated-filter-tab-d', 'disabled'),
         Output('integrated-filter-tab-d', 'style', allow_duplicate=True)],
        [Input('copy2gwas-selected-samples-store', 'data')],
        prevent_initial_call=True
    )
    def disable_filter_d_tab(selected_samples):
        """샘플 수에 따라 Filter D 탭 비활성화"""
        sample_count = len(selected_samples) if selected_samples else 0
        
        if sample_count < 2:
            # 2개 미만 샘플: 탭 비활성화
            disabled = True
            style = {
                'opacity': '0.5',
                'pointer-events': 'none',
                'cursor': 'not-allowed'
            }
        else:
            # 2개 이상 샘플: 탭 활성화
            disabled = False
            style = {
                'opacity': '1.0',
                'pointer-events': 'auto',
                'cursor': 'pointer'
            }
        
        return disabled, style
        '''


setup_integrated_filter_callbacks()


# 🎯 Toggle 및 Searchable Varieties 기능 제거됨 (좌우 분할 레이아웃으로 변경)

# 🎯 노드 클릭 처리 (Expand/Selection 모드) - 선택 품종 유지

def get_selected_stylesheet(selected_nodes):
    """선택된 노드들을 포함한 스타일시트 반환"""
    stylesheet = get_default_stylesheet()
    
    for selected_id in selected_nodes:
        stylesheet.append({
            'selector': f'node[id="{selected_id}"]',
            'style': {
                'background-color': '#ff69b4',
                'border-width': 3,
                'border-color': '#ff1493'
            }
        })
    
    return stylesheet

# ------------------- HELPER -------------------



def get_available_nodes_from_elements(elements, *, need_any=True, want_vcf=True, want_pheno=True):
    """
    Cytoscape elements에서 group='nodes'만 보고
    has_vcf / has_pheno 기준으로 '가용(Available) 노드' 목록을 만든다.
    반환: [{id,name,it_number,has_vcf,vcf_id,has_pheno}, ...]
    """
    print(f"🔍 Available nodes 스캔 시작: {len(elements) if elements else 0}개 elements")
    out, seen = [], set()
    for el in (elements or []):
        if el.get('group', 'nodes') != 'nodes':
            continue
        d = el.get('data', {}) or {}
        nid = d.get('id')
        if not nid or nid in seen:
            continue
        has_v = bool(d.get('has_vcf'))
        # has_pheno는 has_pheno 속성이 True이거나 it_number가 존재하면 True
        itn = d.get('it_number')
        has_p = bool(itn) and itn != "No information"
        
        print(f"  노드 스캔: {nid}, has_vcf={has_v}, it_number={itn}, has_pheno={has_p}")
        
        ok = ((want_vcf and has_v) or (want_pheno and has_p)) if need_any \
             else ((not want_vcf or has_v) and (not want_pheno or has_p))
        if ok:
            seen.add(nid)
            node_data = {
                'id': nid,
                'name': d.get('label') or d.get('name'),
                'it_number': d.get('it_number'),
                'has_vcf': has_v,
                'vcf_id': d.get('vcf_status') if has_v else None,
                'has_pheno': has_p
            }
            out.append(node_data)
            print(f"    ✅ Available 노드 추가: {node_data}")
    
    print(f"🎯 Available nodes 스캔 완료: {len(out)}개 노드")
    return out


def _esc_attr_value(s: str) -> str:
    """Cytoscape CSS 속성 선택자에서 안전하게 쓰기 위한 이스케이프"""
    s = str(s)
    return s.replace("\\", "\\\\").replace('"', '\\"')


def build_stylesheet_with_available(*, base=None, available_ids=None, expanded_ids=None):
    """
    기본 스타일 + '가용 노드(밝은 청록 #5dade2)'와 '확장 노드 링(녹색 border)' 적용.
    Fixed: 일관된 색상 체계와 텍스트 표시 개선 (v2024)
    ※ id에 공백/특수문자가 있어도 동작하도록 속성 선택자 사용:
       selector: f'node[id = "{...}"]'
    """
    base = base or []
    available_ids = list(set(available_ids or []))
    expanded_ids  = list(set(expanded_ids  or []))

    # 기본 스타일시트 가져오기 (일관성을 위해)
    ss = get_default_stylesheet()

    # data-available 색 적용 (밝은 청록색) - 기본 블루와 구분, 가독성 개선
    ss += [
        {'selector': f'node[id = "{_esc_attr_value(nid)}"]',
         'style': {
             'background-color': '#5dade2',  # 밝은 블루 (가독성 개선)
             'border-width': 2,
             'border-color': '#2980b9',
             'color': '#000',  # 텍스트 대비 강화
             'font-size': 12,  # 폰트 크기 약간 축소
             'text-max-width': 65,  # 텍스트 너비 약간 증가
         }}
        for nid in available_ids
    ]

    # expanded 노드에 선명한 링 적용 (최고 우선순위)
    ss += [
        {'selector': f'node[id = "{_esc_attr_value(nid)}"]',
         'style': {
             'border-width': 4,  # 두께 증가
             'border-color': '#27ae60',  # 선명한 녹색
             'border-opacity': 1
         }}
        for nid in expanded_ids
    ]

    return ss + base

def get_default_stylesheet():
    """🚩 개선된 기본 스타일시트 - expand history 지원, 텍스트 표시 개선 (v2024)"""
    return [
        {'selector': 'node', 'style': {
            'content': 'data(label)',
            'text-valign': 'center',
            'text-halign': 'center',
            'background-color': '#6FB1FC',
            'width': 70,  # 녹이 증가로 텍스트 공간 확보
            'height': 70,  # 높이 증가로 텍스트 공간 확보
            'font-size': 11,  # 폰트 크기 약간 축소
            'color': '#000',
            'text-wrap': 'wrap',
            'text-max-width': 68,  # 영역 내에서 텍스트 래핑
            'text-overflow-wrap': 'anywhere',  # 긴 단어 줄바꿈
            'border-width': 2,
            'border-color': '#34495e'
        }},
        # 🎯 검색 결과 노드 (최우선 - pink fill + red stroke)
        {'selector': 'node.search-result', 'style': {
            'background-color': '#ff69b4',
            'border-width': 5,
            'border-color': '#e74c3c'
        }},
        # 🎯 선택된 노드 (.selected.color-mapped와 동일)
        {'selector': 'node.selected-node', 'style': {
            'background-color': '#ff69b4',
            'border-width': 3,
            'border-color': '#ff1493'
        }},
        # 🚩 expand 클릭된 노드 (green stroke)
        {'selector': 'node.expanded', 'style': {
            'border-width': 4,
            'border-color': '#28a745'
        }},
        # 🚩 expand로 추가된 자식 노드 (orange fill)
        {'selector': 'node.trait', 'style': {
            'background-color': '#ffa500',
            'border-width': 2,
            'border-color': '#ff8c00'
        }},
        # 🎯 GWAS 선택된 노드 (dark gray)
        {'selector': 'node.gwas-selected', 'style': {
            'background-color': '#2c3e50',
            'border-width': 3,
            'border-color': '#34495e',
            'color': '#ffffff'
        }},
        # 🎯 pedigree에서 선택된 노드 (green)
        {'selector': 'node.pedigree-selected', 'style': {
            'background-color': '#27ae60',
            'border-width': 3,
            'border-color': '#229954',
            'color': '#ffffff'
        }},
        # 기본 선택 스타일
        {'selector': 'node:selected', 'style': {
            'background-color': '#B7E1A1',
            'border-width': 2,
            'border-color': '#698B69'
        }},
        {'selector': 'edge', 'style': {
            'curve-style': 'bezier',
            'target-arrow-shape': 'triangle',
            'arrow-scale': 1.5,
            'line-color': '#ccc',
            'target-arrow-color': '#ccc'
        }},
        {'selector': 'node:active', 'style': {
            'overlay-opacity': 0,
            'background-color': '#61bffc',
        }}
    ]

# =============================================================================
# MAIN EXECUTION
# =============================================================================
if __name__ == '__main__':
    # 개발/테스트용 독립 실행
    app = dash.Dash(__name__, 
                   external_stylesheets=[dbc.themes.BOOTSTRAP],
                   suppress_callback_exceptions=True)
    
    # 중복 콜백 허용 설정 - allow_duplicate=True with initial calls
    app.config.prevent_initial_callbacks = 'initial_duplicate'
    
    # Register mannequin copy callbacks with this app instance
    #from mannequin_copy_callbacks import register_mannequin_copy_callbacks
    #register_mannequin_copy_callbacks(app)
    
    app.layout = create_layout()  # 👈 파라미터 없이 기본 레이아웃
    app.run_server(debug=True, port=9099)

def get_available_traits(phenotype_data):
    """표현형 데이터에서 사용 가능한 트레이트 목록 추출"""
    if not phenotype_data or 'data' not in phenotype_data:
        return []
    
    df = pd.DataFrame(phenotype_data['data'])
    if df.empty:
        return []
    
    # 'id' 컬럼 제외하고 데이터가 있는 컬럼들만 반환
    traits = []
    for col in df.columns:
        if col != 'id' and not df[col].isna().all():
            traits.append(col)
    
    return traits

def calculate_trait_color_scale(phenotype_data, trait_name):
    """🎯 DEPRECATED: 기존 함수는 Global Scale을 사용하지 않으므로 Enhanced 버전으로 리다이렉트"""
    print(f"🔄 calculate_trait_color_scale → calculate_enhanced_trait_color_scale 리다이렉트")
    
    # 카테고리 자동 감지
    category_name = detect_trait_category(trait_name)
    
    # 🎯 Enhanced 버전 사용 (Global Scale 보장)
    return calculate_enhanced_trait_color_scale(phenotype_data, trait_name, category_name)

def apply_trait_colors_to_stylesheet(base_stylesheet, color_mapping):
    """🎯 개선된 트레이트 컬러 적용 - 새로운 색상 매핑 구조 지원"""
    if not color_mapping:
        print("   ❌ apply_trait_colors_to_stylesheet: color_mapping이 비어있음")
        return base_stylesheet
    
    print(f"   🎨 apply_trait_colors_to_stylesheet: {len(color_mapping)}개 품종에 색상 적용")
    enhanced_stylesheet = base_stylesheet.copy()
    
    for variety_id, color_info in color_mapping.items():
        # 🎯 새로운 구조: color_info가 dict인지 string인지 확인
        if isinstance(color_info, dict):
            color = color_info['color']
            trait_value = color_info.get('value', 'N/A')
            trait_type = color_info.get('trait_type', 'unknown')
            print(f"      🌾 {variety_id} → {color} (값: {trait_value}, 타입: {trait_type})")
        else:
            # 기존 호환성 유지 (단순 색상 문자열)
            color = color_info
            print(f"      🌾 {variety_id} → {color} (기존 형식)")
        
        # 기존 스타일에서 해당 노드 스타일 제거 (중복 방지)
        enhanced_stylesheet = [style for style in enhanced_stylesheet 
                             if style.get('selector') != f'node[id="{variety_id}"]']
        
        # 🎯 Trait 값에 따른 전체 노드 색상 적용 (배경색을 trait 색상으로 변경)
        enhanced_stylesheet.append({
            'selector': f'node[id="{variety_id}"]',
            'style': {
                # 🌈 배경색을 trait 값에 따른 색상으로 설정 (명확한 시각화)
                'background-color': color,      # trait 값에 따른 색상 적용
                'background-opacity': 0.9,      # 약간의 투명도로 가독성 향상
                
                # 테두리는 더 진한 색상으로 강조
                'border-width': '3px',          
                'border-color': '#ffffff',      # 흰색 테두리로 노드 구분
                'border-style': 'solid',
                'border-opacity': 1,
                
                # 노드 크기 및 z-index
                'width': '75px',
                'height': '75px',
                'z-index': 999,
                
                # 텍스트 가독성 향상 (색상에 따라 적응적 텍스트 색상)
                'color': '#ffffff',
                'text-outline-width': '2px',
                'text-outline-color': '#000000',
                'font-weight': 'bold',
                'font-size': '12px',
                
                # 선택 상태에서도 스타일 유지
                'overlay-opacity': 0,
                'overlay-color': 'transparent',
                
                # 그림자 효과로 입체감 부여
                'box-shadow': '0 4px 8px rgba(0,0,0,0.3)'
            }
        })
    
    print(f"   📊 최종 스타일시트 크기: {len(enhanced_stylesheet)}개 규칙")
    return enhanced_stylesheet


# =============================================================================
# CATEGORY MAPPINGS WITH 6-ITEM GROUPING
# =============================================================================
def get_all_category_mappings():
    """카테고리 분할 없이 원본 매핑을 그대로 반환"""
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

# 기존 호환성을 위한 원본 카테고리 매핑
CATEGORY_MAPPINGS = {
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

# =============================================================================
# PHENOTYPE PLOT UPDATE CALLBACKS
# =============================================================================
'''
@callback(
    [Output(f'phenotype-plot-{category}', 'children') for category in get_all_category_mappings().keys()] +
    [Output(f'trait-spans-{category}', 'children') for category in get_all_category_mappings().keys()],
    [Input('phenotype-data-store', 'data'),
     #Input('search-selected-variety-name', 'data')
     ],
    [State('fixed-vcf-item', 'data'),
     State('fixed-vcf-item-nopedi', 'data')],
    prevent_initial_call=True
)
def update_all_phenotype_plots(applied_nodes,  fixed_vcf_value, fixed_vcf_value_nopedi): #search_variety_name,
    """선택된 품종이 변경될 때 모든 표현형 플롯과 trait spans 업데이트 — DB 제거 버전"""

    #applied_data = applied_nodes or []
    #selected_nodes = [node.get('it_number') for node in applied_data if node.get('it_number')] if applied_data else []
    applied_data = applied_nodes or []
    selected_nodes = [node.get('id') for node in applied_data if node.get('id')]  # ✅ 수정된 부분

    nopedi_variety_id = (fixed_vcf_value_nopedi or {}).get('variety_id')
    pedi_variety_id   = (fixed_vcf_value or {}).get('variety_id')

    # 분석 대상 노드 결정
    #print(f"search_variety_name: {search_variety_name}")
    #print(f"nopedi_variety_id: {nopedi_variety_id}")
    #print(f"pedi_variety_id: {pedi_variety_id}")
    #print(f"selected_nodes: {selected_nodes}")
    #if nopedi_variety_id and search_variety_name:
    #    analysis_nodes = list(dict.fromkeys([nopedi_variety_id, search_variety_name]))
    if selected_nodes:
        analysis_nodes = selected_nodes
    elif pedi_variety_id:
        analysis_nodes = [pedi_variety_id]
    else:
        # 아무 것도 없으면 빈 컨테이너 반환
        cats = list(get_all_category_mappings().keys())
        empty_containers = [html.Div([html.P("No varieties selected",
                              style={'textAlign': 'center', 'color': '#666', 'margin': '50px'})])
                            for _ in cats]
        empty_spans = [html.Div() for _ in cats]
        return empty_containers + empty_spans

    # 선택된 품종들의 표현형 데이터 조회 (PHENO_DF 기반)
    #print(f"analysis_nodes: {analysis_nodes}")
    selected_varieties_data = {}
    for node_id in analysis_nodes:
        try:
            # 이 단계에서 node_id == IT Number 로 통일(이미 그렇게 쓰고 계신 흐름 유지)
            it_number = node_id
            if not it_number or it_number == 'N/A':
                continue

            df = pheno_by_id(it_number)
            if not df.empty:
                selected_varieties_data[node_id] = df.iloc[0].to_dict()
        except Exception:
            continue
    #print(f"df.head(): {df.head()}")
    #print(f"selected_varieties_data: {selected_varieties_data}")

    all_category_mappings = get_all_category_mappings()

    if not selected_varieties_data:
        empty_containers = [html.Div([html.P("No phenotype data found for selected varieties",
                              style={'textAlign': 'center', 'color': '#666', 'margin': '50px'})])
                            for _ in all_category_mappings.keys()]
        empty_spans = [html.Div() for _ in all_category_mappings.keys()]
        return empty_containers + empty_spans

    # 각 카테고리별 플롯/스팬 생성
    containers = []
    spans = []
    for category_name in all_category_mappings.keys():
        try:
            plot_container, trait_spans = _category_plots_from_pheno_df(category_name, selected_varieties_data, selected_varieties_data or {})
            containers.append(plot_container)
            # 별도 spans 영역을 사용하신다면 trait_spans 배치
            # (기존 코드에서는 빈 컨테이너로 대체했으므로, 필요 시 아래 라인 교체)
            spans.append(html.Div(trait_spans))
        except Exception as e:
            containers.append(html.Div([
                html.P(f"Error creating plot for {category_name}: {str(e)}",
                       style={'textAlign': 'center', 'color': 'red', 'margin': '50px'})
            ]))
            spans.append(html.Div())

    return containers + spans
'''

# =============================================================================
# TRAIT COLOR SCHEME CALLBACKS
# =============================================================================

# 기존 드롭다운 콜백들 제거됨 - 버튼 방식으로 대체

# 트레이트 선택 버튼 생성 콜백들 (2개 이상 선택 시에만 활성화)
for category in CATEGORY_MAPPINGS.keys():
    @callback(
        Output(f'trait-selection-buttons-{category}', 'children'),
        [Input('phenotype-data-store', 'data')],  # Applied 데이터에 반응 (unified for both regular and nopedi cases)
        [State('fixed-vcf-item', 'data')],  # fixed variety를 State로 추가
        prevent_initial_call=True  # 🎯 초기 로딩도 허용하여 첫 진입 시에도 버튼 표시
    )
    def update_trait_selection_buttons(applied_nodes, fixed_vcf_value, category_name=category):
        """카테고리별 트레이트 선택 버튼 생성 (2개 이상 선택된 경우에만)"""
        # applied_nodes - now unified for both regular and nopedi cases
        applied_data = applied_nodes or []
        # applied_data는 [{id, name, ...}, ...] 구조이므로 ID만 추출
        selected_nodes = [node.get('id') for node in applied_data if node.get('id')] if applied_data else []
        
        # fixed-variety-item이 있고 applied가 없으면 추가
        if not selected_nodes and fixed_vcf_value and fixed_vcf_value != '-':
            selected_nodes = [fixed_vcf_value]
            
        if not selected_nodes or len(selected_nodes) < 2:
            return html.Div()  # 빈 컨테이너 반환
        
        try:
            # IT번호가 있는 품종들의 표현형 데이터만 수집
            pedigree_app = get_pedigree_app()
            phenotype_data = {'data': []}
            it_nodes_count = 0
            
            for node_id in selected_nodes:
                try:
                    variety_info = get_variety_info_by_pedigree_cached(str(node_id))
                    it_number = variety_info.get("IT Number")
                    
                    if it_number and it_number != 'N/A':
                        it_nodes_count += 1
                        df = get_db_data(it_number)
                        if not df.empty:
                            row_data = df.iloc[0].to_dict()
                            row_data['id'] = node_id
                            phenotype_data['data'].append(row_data)
                except Exception as e:
                    continue
            
            if not phenotype_data['data']:
                return html.Div()  # 빈 컨테이너 반환 (Trait Mapping Available 제거)
            
            # 해당 카테고리의 트레이트들만 필터링
            all_category_mappings = get_all_category_mappings()
            category_traits = all_category_mappings[category_name]['traits']
            available_traits = get_available_traits(phenotype_data)
            valid_traits = [trait for trait in category_traits if trait in available_traits]
            
            if not valid_traits:
                return html.Div()  # 빈 컨테이너 반환 (Trait Mapping Available 제거)
            
            # 트레이트 버튼들 생성
            buttons = []
            trait_icons = {
                'numeric': '📊',
                'categorical': '🏷️',
                'date': '📅'
            }
            
            for trait in valid_traits:
                trait_type = get_trait_type_from_db(trait)
                icon = trait_icons.get(trait_type, '📊')
                
                button = dbc.Button(
                    f"{icon} {trait.replace('_', ' ').title()}",
                    id={'type': 'trait-color-button', 'trait': trait, 'category': category_name},
                    color="outline-primary",
                    size="sm",
                    className="me-2 mb-2"
                )
                buttons.append(button)
            
            return html.Div()  # Trait Mapping Available 영역 완전 제거
            
        except Exception as e:
            print(f"Error creating trait selection buttons for {category_name}: {e}")
            return html.Div(f"Error: {str(e)}", className="text-danger")

# 트레이트 span과 버튼 클릭 콜백 통합
@callback(
    [Output('pedigree-cytoscape', 'stylesheet', allow_duplicate=True),
     Output('trait-legend-container', 'children'),
     Output('trait-legend-section', 'style'),
     Output('trait-color-store', 'data'),
    Output('active-trait-store', 'data')  # ✅ 추가: 현재 클릭 상태 저장
      ] +
    [Output(f'phenotype-plot-{category}', 'children', allow_duplicate=True)
     for category in get_all_category_mappings().keys()],
    [Input({'type': 'trait-title-button', 'trait': ALL, 'category': ALL}, 'n_clicks'),
     Input({'type': 'trait-color-button', 'trait': ALL, 'category': ALL}, 'n_clicks'),
     Input({'type': 'trait-name-span',   'trait': ALL, 'category': ALL}, 'n_clicks')],
    [State('pedigree-cytoscape', 'elements'),   # ✅ 전체 노드/엣지만 사용
     State('trait-color-store', 'data'),
     State('child-nodes-store', 'data'),
     State('highlighted-components-store', 'data'),
    State('component-info-store', 'data'),
    State('active-trait-store', 'data'),
    State('selected-nodes-store', 'data'),
    ],  # ✅ 이전 클릭 정보
    
    prevent_initial_call=True
)
def apply_trait_color_from_interactions(
    title_button_clicks, button_clicks, span_clicks,
    elements,                    # ✅ selectedNodeData 제거
    current_colors, child_nodes,selected_comps, comp_info,#, active_trait
    active_trait_store,selected_nodes                                                    
):
    print("🚩 Callback triggered")

    ctx = dash.callback_context
    # 버튼 중 하나라도 클릭된 게 없으면 업데이트 중단
    #print(ctx.inputs_list)
    if not ctx.triggered or (not any(title_button_clicks or []) and
                             not any(button_clicks or []) and
                             not any(span_clicks or [])):
        raise PreventUpdate

    # --- 어떤 입력이 트리거됐는지 파싱 ---
    clicked_trait, clicked_category = None, None
    try:
        if title_button_clicks:
            for i, c in enumerate(title_button_clicks):
                if c:
                    iid = ctx.inputs_list[0][i]['id']
                    clicked_trait   = iid['trait']
                    clicked_category = iid['category']
                    break
        if not clicked_trait and button_clicks:
            for i, c in enumerate(button_clicks):
                if c:
                    iid = ctx.inputs_list[1][i]['id']
                    clicked_trait   = iid['trait']
                    clicked_category = iid['category']
                    break
        if not clicked_trait and span_clicks:
            for i, c in enumerate(span_clicks):
                if c:
                    iid = ctx.inputs_list[2][i]['id']
                    clicked_trait   = iid['trait']
                    clicked_category = iid['category']
                    break
    except Exception as e:
        print(f"❌ Error parsing trigger: {e}")

    if not clicked_trait:
        print("🚫 No trait parsed from trigger")
        raise PreventUpdate

    previous_trait = active_trait_store.get("trait") if active_trait_store else None
    previous_category = active_trait_store.get("category") if active_trait_store else None
    # --- 같은 버튼 다시 클릭 → OFF 처리 ---
    if previous_trait == clicked_trait and previous_category == clicked_category:
        print(f"🔁 Same trait clicked again → OFF mode")
        new_active = {"trait": None, "category": None}

        # --- elements 에서 모든 노드의 id/IT Number 수집 ---
        node_ids, it_numbers, node_to_itn = _extract_node_it_map(elements)
        print(f"📌 extracted nodes: {len(node_ids)} / ITs: {len(it_numbers)}")
        if not node_to_itn:
            print("🚫 No IT numbers found in elements")
            raise PreventUpdate

        # --- PHENO_DF subset ---
        df_subset = pheno_by_ids(it_numbers)   # (id 열이 IT Number인 DF 반환 가정)
        if df_subset.empty or clicked_trait not in df_subset.columns:
            print("🚫 df_subset empty or trait not in columns")
            raise PreventUpdate

        # phenotype_data 구성 (node id로 재매핑)
        phenotype_data = {'data': []}
        selected_varieties_data = {}
        # df_subset에서 id(=IT Number) → 행 매핑
        by_itn = {str(r['id']): r for _, r in df_subset.iterrows() if 'id' in r}

        for nid, itn in node_to_itn.items():
            row = by_itn.get(str(itn))
            if row is not None:
                row_dict = dict(row)
                row_dict['id'] = nid  # cytoscape node id로 치환
                phenotype_data['data'].append(row_dict)
                selected_varieties_data[nid] = row

        if not phenotype_data['data']:
            print("🚫 phenotype_data empty after remap")
            raise PreventUpdate

        # --- 색상 매핑 ---
        detected_category = detect_trait_category(clicked_trait)
        color_mapping = calculate_enhanced_trait_color_scale(
            phenotype_data, clicked_trait, detected_category
        )
        if not color_mapping:
            print("🚫 color_mapping empty")
            raise PreventUpdate

        # Cytoscape 초기화
        base_stylesheet = get_default_stylesheet2()
        highlight_style = []
        if selected_comps:
            color_palette = [
                "#8e44ad", "#6c5ce7", "#2d3436", "#1a1a1a",
                "#c0399f", "#be90d4", "#4b0082", "#5c3c92", "#3d3d3d",
                "#1C6BA0", "#2980b9", "#e74c3c", "#c0392b", "#e67e22",
                "#d35400", "#27ae60", "#2ecc71", "#f1c40f", "#f39c12"
            ]
            for i, comp_id in enumerate(selected_comps):
                comp = next((c for c in comp_info if c.get('id') == comp_id), None)
                if not comp:
                    continue
                ccolor = color_palette[i % len(color_palette)]
                
                for e in comp.get("edges", []):
                    src, tgt = e.get("source"), e.get("target")
                    if src and tgt:
                        highlight_style.append({
                            "selector": f'edge[source="{src}"][target="{tgt}"]',
                            "style": {
                                "line-color": ccolor,
                                "target-arrow-color": ccolor,
                                "width": 4,
                            },
                        })
        highlight_select_styles = []
        for nd in selected_nodes or []:
            nid = nd.get("id") or nd.get("name")
            if not nid:
                continue
            highlight_select_styles.append({
                "selector": f'node[id = "{nid}"]',
                "style": {
                    "border-color": "#ff33cc",   # 짙은 초록 외곽선
                    "border-width": 6,           # 두꺼운 강조선
                    "transition-property": "border-color, border-width",
                    "transition-duration": "0.3s",
                },
            })

        reset_stylesheet = base_stylesheet + highlight_style+highlight_select_styles

        # legend 숨기기
        legend = create_enhanced_trait_legend(
            phenotype_data, clicked_trait, color_mapping, detect_trait_category(clicked_trait)
        )
        legend_style = {"display": "none"}
        color_mapping = {}

        # ✅ 여기서도 plot regeneration 그대로 수행
        updated_plot_containers = []
        for category_name in get_all_category_mappings().keys():
            try:
                plot_container = _category_plots_from_pheno_df(
                    category_name, PHENO_DF, selected_varieties_data,
                    new_active["trait"], new_active["category"]   # None, None 전달
                )
                print(f"✅ Plot updated (OFF): {category_name}")
            except Exception as e:
                print(f"❌ Plot error (OFF) {category_name}: {e}")
                plot_container = html.Div()
            updated_plot_containers.append(plot_container)

        # 반환 (OFF 상태)
        return [reset_stylesheet, legend, legend_style, color_mapping, new_active] + updated_plot_containers
        # --- 버튼 토글 상태 결정 ---
    
    new_active = {"trait": clicked_trait, "category": clicked_category}
    

    

    # --- elements 에서 모든 노드의 id/IT Number 수집 ---
    node_ids, it_numbers, node_to_itn = _extract_node_it_map(elements)
    print(f"📌 extracted nodes: {len(node_ids)} / ITs: {len(it_numbers)}")
    if not node_to_itn:
        print("🚫 No IT numbers found in elements")
        raise PreventUpdate

    # --- PHENO_DF subset ---
    df_subset = pheno_by_ids(it_numbers)   # (id 열이 IT Number인 DF 반환 가정)
    if df_subset.empty or clicked_trait not in df_subset.columns:
        print("🚫 df_subset empty or trait not in columns")
        raise PreventUpdate

    # phenotype_data 구성 (node id로 재매핑)
    phenotype_data = {'data': []}
    selected_varieties_data = {}
    # df_subset에서 id(=IT Number) → 행 매핑
    by_itn = {str(r['id']): r for _, r in df_subset.iterrows() if 'id' in r}

    for nid, itn in node_to_itn.items():
        row = by_itn.get(str(itn))
        if row is not None:
            row_dict = dict(row)
            row_dict['id'] = nid  # cytoscape node id로 치환
            phenotype_data['data'].append(row_dict)
            selected_varieties_data[nid] = row

    if not phenotype_data['data']:
        print("🚫 phenotype_data empty after remap")
        raise PreventUpdate

    # --- 색상 매핑 ---
    detected_category = detect_trait_category(clicked_trait)
    color_mapping = calculate_enhanced_trait_color_scale(
        phenotype_data, clicked_trait, detected_category
    )
    if not color_mapping:
        print("🚫 color_mapping empty")
        raise PreventUpdate

    # --- Cytoscape stylesheet 생성 ---
    base_stylesheet = get_default_stylesheet2()#get_enhanced_stylesheet(node_ids, [], child_nodes or [])
    trait_specific_styles = []
    for node_id in node_ids:
        cinfo = color_mapping.get(str(node_id))
        trait_color = (cinfo.get('color') if isinstance(cinfo, dict) else cinfo) or '#bdc3c7'
        if trait_color:
            trait_specific_styles.append({
                'selector': f'node[id="{node_id}"]',
                'style': {
                    'background-color': trait_color,
                    'border-color': trait_color,
                    'border-width': '3px'
                }
            })

    highlight_style = []
    if selected_comps:
        color_palette = [
            "#8e44ad",  "#6c5ce7", "#2d3436", "#1a1a1a",
            "#c0399f", "#be90d4", "#4b0082", "#5c3c92", "#3d3d3d",
            "#1C6BA0", "#2980b9", "#e74c3c", "#c0392b", "#e67e22",
            "#d35400", "#27ae60", "#2ecc71", "#f1c40f", "#f39c12"
        ]
        for i, comp_id in enumerate(selected_comps):
            comp = next((c for c in comp_info if c.get('id') == comp_id), None)
            if not comp:
                continue
            ccolor = color_palette[i % len(color_palette)]
            for e in comp.get("edges", []):
                src, tgt = e.get("source"), e.get("target")
                if src and tgt:
                    highlight_style.append({
                        "selector": f'edge[source="{src}"][target="{tgt}"]',
                        "style": {
                            "line-color": ccolor,
                            "target-arrow-color": ccolor,
                            "width": 4,
                        },
                    })
    highlight_select_styles = []
    for nd in selected_nodes or []:
        nid = nd.get("id") or nd.get("name")
        if not nid:
            continue
        highlight_select_styles.append({
            "selector": f'node[id = "{nid}"]',
            "style": {
                "border-color": "#ff33cc",   # 짙은 초록 외곽선
                "border-width": 6,           # 두꺼운 강조선
                "transition-property": "border-color, border-width",
                "transition-duration": "0.3s",
            },
        })

    enhanced_stylesheet = base_stylesheet + trait_specific_styles + highlight_style+highlight_select_styles

    # --- 범례 ---
    legend = create_enhanced_trait_legend(
        phenotype_data, clicked_trait, color_mapping, detected_category#clicked_category
    )
    legend_style = {'padding': '2px', 'display': 'block',
                                            "width": "100%",                   # ✅ 전체 폭 비율 조절
                                            "marginLeft": "auto",             # ✅ 오른쪽 정렬
                                            "marginRight": "0",
                                            "backgroundColor": "#f8f9fa",
                                            "border": "1px solid #dee2e6",
                                            "borderRadius": "8px",
                                            "padding": "10px 14px",
                                            "boxShadow": "0 2px 6px rgba(0,0,0,0.1)",
                                            "transition": "width 0.3s ease-in-out",
                                            "marginTop": "8px",
                                            "marginBottom": "6px",}

    # --- 카테고리별 플롯 업데이트 ---
    updated_plot_containers = []
    for category_name in get_all_category_mappings().keys():
        try:
            plot_container = _category_plots_from_pheno_df(
                category_name, PHENO_DF, selected_varieties_data,new_active["trait"], new_active["category"]
            )
            print(f"✅ Plot updated: {category_name}")
        except Exception as e:
            print(f"❌ Plot error {category_name}: {e}")
            plot_container = html.Div()
        updated_plot_containers.append(plot_container)

    print("🎉 Callback success")
    return [enhanced_stylesheet, legend, legend_style, color_mapping, new_active] + updated_plot_containers

def create_enhanced_trait_legend(phenotype_data, trait_name, color_mapping, category_name):
    """PHENO_DF 기반 범례 생성: numeric, categorical/date 구분 (마커/품종 상세 포함)"""
    import pandas as pd
    from dash import html

    print("\n🚩 === create_enhanced_trait_legend CALLED ===")
    print(f"🔹 trait_name = {trait_name}, category_name = {category_name}")
    print(f"🔹 phenotype_data keys = {list(phenotype_data.keys()) if phenotype_data else None}")
    print(f"🔹 color_mapping size = {len(color_mapping) if color_mapping else 0}")

    if not phenotype_data or "data" not in phenotype_data:
        return html.Div()

    df = pd.DataFrame(phenotype_data["data"])
    if df.empty or trait_name not in df.columns:
        return html.Div()

    theme = get_category_color_theme(category_name)
    trait_type = determine_trait_type_fallback(trait_name)
    legend_items = []

    # ---------------- Numeric Case ----------------
    if trait_type == "numeric":
        global_min = global_max = None
        global_count = global_avg = 0
        selected_variety_info = []

        for vid, cinfo in color_mapping.items():
            if isinstance(cinfo, dict):
                if global_min is None:
                    global_min = cinfo.get("global_min")
                    global_max = cinfo.get("global_max")
                    global_count = cinfo.get("global_count", 0)
                    global_avg = cinfo.get("global_avg", 0)
                selected_variety_info.append({
                    "id": vid,
                    "value": cinfo.get("value", 0),
                    "color": cinfo.get("color", "#666"),
                    "normalized": cinfo.get("normalized", 0.5),
                })

        # fallback global range
        if global_min is None or global_max is None:
            vals = pd.to_numeric(df[trait_name], errors="coerce").dropna()
            if not vals.empty:
                global_min, global_max = float(vals.min()), float(vals.max())
                global_count, global_avg = int(vals.size), float(vals.mean())
            else:
                global_min, global_max, global_count, global_avg = 0.0, 100.0, 0, 50.0

        # Title
        legend_items.append(html.Div(f"🎨 {trait_name}",
            style={"fontWeight":"bold","fontSize":"14px","marginBottom":"8px","color":theme["primary"]}))

        # Gradient bar
        def _css_gradient(steps=10):
            colors = [generate_enhanced_gradient_color(i/steps, theme, category_name) for i in range(steps+1)]
            return f"linear-gradient(to right, {', '.join(colors)})"

        gradient_style = _css_gradient(10)
        min_label = html.Div(f"{global_min:.2f}", style={"fontSize":"12px","color":"#666","minWidth":"40px","textAlign":"right"})
        max_label = html.Div(f"{global_max:.2f}", style={"fontSize":"12px","color":"#666","minWidth":"40px","textAlign":"left"})
        gradient_bar = html.Div(style={"width":"100%","height":"10px","background":gradient_style,
                                       "border":"1px solid #ddd","borderRadius":"4px"})

        # Marker overlay
        marker_items = []
        for variety in selected_variety_info:
            value = variety.get("value")
            norm = variety.get("normalized")
            if value is None:
                continue
            if global_max > global_min:
                nv = norm if norm is not None else (value - global_min) / (global_max - global_min)
                nv = max(0.0, min(1.0, nv))
                left_px = nv * 200
            else:
                left_px = 100

            marker = html.Div(
                [
                    # Vertical line
                    html.Div(style={"width":"2px","height":"20px","backgroundColor":variety["color"],
                                    "border":"1px solid white","borderRadius":"1px"}),
                    # Top dot
                    html.Div(style={"width":"6px","height":"6px","backgroundColor":variety["color"],
                                    "border":"1px solid white","borderRadius":"50%","marginTop":"-3px","marginLeft":"-2px"}),
                    # Tooltip
                    html.Div(f"{variety['id']}: {value:.2f}",
                        style={"position":"absolute","bottom":"25px","left":"50%","transform":"translateX(-50%)",
                               "backgroundColor":"rgba(0,0,0,0.8)","color":"white","padding":"2px 6px",
                               "borderRadius":"4px","fontSize":"10px","opacity":"0","pointerEvents":"none"})
                ],
                style={"position":"absolute","left":f"{left_px}px","top":"0px",
                       "display":"flex","flexDirection":"column","alignItems":"center","cursor":"pointer"},
                className="gradient-marker"
            )
            marker_items.append(marker)

        gradient_container = html.Div(
            [gradient_bar, html.Div(marker_items, style={"position":"absolute","top":"0px","left":"0px",
                                                         "width": "100%","height":"10px"})],
            style={"position":"relative","width":"100%","height":"10px"}
        )

        legend_items.append(
            html.Div([min_label, gradient_container, max_label],
                     style={"display":"flex","alignItems":"center","marginBottom":"12px","gap":"8px"})
        )

        # 선택된 품종 리스트
        if selected_variety_info:
            legend_items.append(html.Hr(style={"margin":"10px 0","opacity":"0.4"}))
            legend_items.append(html.Div("📍 Variety with phenotype data",
                                         style={"fontSize":"12px","fontWeight":"bold","color":theme["primary"]}))
            for variety in selected_variety_info:
                legend_items.append(
                    html.Div([
                        html.Div("●", style={"color":variety["color"],"fontSize":"18px","marginRight":"6px"}),
                        html.Span(f"{variety['id']} – {variety['value']:.2f}", style={"fontSize":"11px","color":"#555"})
                    ], style={"display":"flex","alignItems":"center","marginBottom":"4px"})
                )

    # ---------------- Categorical/Date Case ----------------
    elif trait_type in ["categorical","date"]:
        legend_items.append(html.Div(f"🏷️ {trait_name}",
            style={"fontWeight":"bold","fontSize":"14px","marginBottom":"8px","color":theme["primary"]}))

        unique_categories = {}
        for vid, cinfo in color_mapping.items():
            if isinstance(cinfo, dict):
                val = cinfo.get("value","N/A")
                col = cinfo.get("color","#666")
                if val not in unique_categories:
                    unique_categories[val] = {"color":col,"varieties":[vid]}
                else:
                    unique_categories[val]["varieties"].append(vid)

        # 카테고리별 legend
        for val, info in unique_categories.items():
            legend_items.append(
                html.Div([
                    html.Div(style={"width":"16px","height":"16px","backgroundColor":info["color"],
                                    "border":"1px solid #ddd","borderRadius":"3px","marginRight":"8px"}),
                    html.Span(f"{val}", style={"fontSize":"12px","color":"#555"})
                ], style={"display":"flex","alignItems":"center","marginBottom":"4px"})
            )
            # 품종 리스트
            for vid in info["varieties"]:
                legend_items.append(
                    html.Div(f" - {vid}", style={"fontSize":"10px","color":"#777","marginLeft":"24px"})
                )

    return html.Div(
        legend_items,
        style={"padding":"12px","border":"1px solid #e9ecef","borderRadius":"8px",
               "backgroundColor":"#f8f9fa","marginTop":"10px"}
    )



# =============================================================================

# Phenotype CSV download callback
'''@callback(
    Output("download-phenotype", "data"),
    Input("download-phenotype-csv", "n_clicks"),
    State('multi-selected-nodes', 'data'),
    prevent_initial_call=True
)
def download_phenotype_csv(n_clicks, selected_nodes):
    """선택된 품종들의 표현형 데이터를 CSV로 다운로드"""
    if not n_clicks or not selected_nodes:
        return no_update
    
    try:
        # 선택된 품종들의 표현형 데이터 조회
        selected_varieties_data = {}
        pedigree_app = get_pedigree_app()
        
        for node_id in selected_nodes:
            # 1단계: 품종명으로부터 IT 번호 가져오기
            try:
                variety_info = get_variety_info_by_pedigree_cached(str(node_id))
                it_number = variety_info.get("IT Number")
                
                if not it_number or it_number == 'N/A':
                    continue
                    
                # 2단계: IT 번호로 표현형 데이터 조회
                df = get_db_data(it_number)
                if not df.empty:
                    selected_varieties_data[node_id] = df.iloc[0].to_dict()
            except Exception as e:
                print(f"Error processing {node_id}: {e}")
                continue
        
        # 데이터프레임으로 변환
        if selected_varieties_data:
            df_export = pd.DataFrame.from_dict(selected_varieties_data, orient='index')
            df_export.index.name = 'variety_id'
            df_export = df_export.reset_index()
            
            # 현재 시간을 파일명에 포함
            import datetime
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"phenotype_data_{timestamp}.csv"
            
            return dcc.send_data_frame(df_export.to_csv, filename, index=False)
        
    except Exception as e:
        print(f"Error in phenotype CSV download: {e}")
        return no_update
    
    return no_update

# 🚀 데이터베이스 쿼리 캐시 - 중복 조회 방지
'''


# ====== Subtrait Color Mapping (고정) ======
ALL_SUBTRAITS = [
    'Yield(14)',  'Stress(66)', 'Plantvigor(23)', 'Biochemical(103)', 'Plantgrowth_Development(17)', 'Plant_quality(26)', 'Biologicalprocess(19)', 'Plantmorphology(191)', 'Sterility_Fertility(2)'
]
ALL_SUBTRAIT_COLORS = px.colors.qualitative.Set1[:9]
SUBTRAIT_COLOR_MAP = {sub: color for sub, color in zip(ALL_SUBTRAITS, ALL_SUBTRAIT_COLORS)}
SUBTRAIT_COLOR_MAP['default'] = '#888888'  # 기본 색상 추가

# =============================================================================
# 🎨 ENHANCED COLOR STRATEGY FUNCTIONS
# =============================================================================


# =============================================================================
# 🔍 PEDIGREE VARIETY SEARCH CALLBACKS (for 'nopedi' cases)
# =============================================================================
'''
# 초기 로딩 시 기본 옵션 설정 콜백
@callback(
    Output('pedigree-variety-search', 'options'),
    Input('pedigree-variety-search', 'id'),  # 컴포넌트 로드 시 트리거
    prevent_initial_call=False               # 초기 로딩 시 옵션 로드 허용
)
def init_pedigree_variety_search_options(_component_id):
    """
    초기 로딩 시 기본 옵션 설정 (첫 20개 표시)
    """
    try:
        print("🔍 Pedigree search 초기 옵션 로드 시작...")
        all_options = get_pedigree_variety_options()
        print(f"✅ Pedigree search 옵션 로드 완료: {len(all_options)}개")
        return all_options[:20]  # 초기에는 첫 20개만 표시
    except Exception as e:
        print(f"❌ 초기 옵션 로드 오류: {e}")
        return []

# 자동완성 옵션 업데이트 콜백 - 값 유지 기능 개선

@callback(
    [Output('pedigree-variety-search', 'options', allow_duplicate=True),
     Output('pedigree-variety-search', 'value', allow_duplicate=True)],
    Input('pedigree-variety-search', 'search_value'),
    State('pedigree-variety-search', 'value'),
    prevent_initial_call=True
)
def update_pedigree_variety_options(search_value, current_value):
    """
    검색어에 따른 계보도 품종 자동완성 옵션 업데이트 - 선택값 유지
    (DB 미사용, PedigreeApp의 CSV 데이터 기반)
    """
    print(f"🔍 자동완성 검색 트리거: search_value='{search_value}', current_value='{current_value}'")

    # 전체 옵션은 매번 PedigreeApp에서 가져오되, CSV는 이미 메모리에 로드된 상태(캐시/싱글톤)여야 합니다.
    all_options = get_pedigree_variety_options()

    # 검색어가 비었으면 기본 20개 + 현재 선택 유지
    if not search_value or len(str(search_value)) < 1:
        print("   검색어 없음 또는 너무 짧음, 기본 옵션 반환")
        default_options = all_options[:20]

        # 현재 선택값 유지 로직
        if current_value:
            if not any(opt['value'] == current_value for opt in default_options):
                current_opt = next((opt for opt in all_options if opt['value'] == current_value), None)
                if current_opt:
                    default_options.insert(0, current_opt)
            return default_options, current_value

        return default_options, None

    # 검색어 필터링 (부분 일치, 최대 50개)
    needle = str(search_value).lower()
    filtered = [opt for opt in all_options if needle in str(opt['label']).lower()][:50]
    print(f"   필터링 결과: {len(filtered)}개")
    if filtered:
        print(f"   첫 번째 결과: {filtered[0]['label']}")

    # 현재 선택값이 필터링 결과에 있으면 유지
    preserved_value = None
    if current_value and any(opt['value'] == current_value for opt in filtered):
        preserved_value = current_value
        print(f"   현재 선택값 유지: {preserved_value}")
    elif current_value:
        print(f"   현재 선택값 필터링됨: {current_value}")

    return filtered, preserved_value

# 검색 버튼 활성화/비활성화 콜백
@callback(
    Output('pedigree-search-button', 'disabled'),
    Input('pedigree-variety-search', 'value'),
    prevent_initial_call=True
)
def update_search_button_state(selected_variety):
    """품종 선택 시 검색 버튼 활성화"""
    print(f"🔍 자동완성 값 선택: {selected_variety}")
    is_disabled = not bool(selected_variety)
    print(f"   버튼 비활성화: {is_disabled}")
    return is_disabled

# ------------------- CALLBACK: build_initial_pedigree -------------------

@callback(
    [Output('pedigree-cytoscape', 'elements'),
     Output('pedigree-cytoscape', 'stylesheet'),
     Output('pedigree-path-store', 'data')],
    Input('pedigree-search-button', 'n_clicks'),
    State('pedigree-variety-search', 'value'),
    prevent_initial_call=True
)
def build_initial_pedigree(_, query):
    """
    - 검색어로 계보 elements 생성(노드 data: id,name/label,it_number,has_vcf,vcf_status,has_pheno,parents,children)
    - 검색 노드에는 classes='search-node' 부여(생성부에서 처리)
    - Expanded 초기화([])
    - stylesheet: 가용 노드(파랑) + expanded 링(초기엔 없음)
    """
    if not query:
        raise exceptions.PreventUpdate
        
    try:
        pedigree_app = get_pedigree_app()
        if pedigree_app:
            # nopedi의 경우 품종명으로 계보도 노드/엣지 생성
            try:
                nodes, edges = pedigree_app.get_connected_nodes(query, 2, 2)  # 부모 2세대, 자식 2세대
                if nodes:
                    elements = pedigree_app.create_cytoscape_elements(nodes, edges)
                    print(f"✅ nopedi: {query} 품종의 계보도 생성 성공 - {len(elements)}개 요소")
                else:
                    print(f"❌ No pedigree data found for '{query}'. Please try another variety.")
                    elements = []
            except Exception as inner_e:
                print(f"❌ No pedigree data found for '{query}'. Please try another variety. Error: {inner_e}")
                elements = []
        else:
            print("❌ PedigreeApp not available")
            elements = []
            
        available = get_available_nodes_from_elements(elements)
        stylesheet = build_stylesheet_with_available(
            base=[],
            available_ids=[x['id'] for x in available],
            expanded_ids=[]
        )
        # pedigree-path-store에 검색된 품종명 추가
        path_data = [query] if query else []
        return elements, stylesheet, path_data
    except Exception as e:
        print(f"Error in build_initial_pedigree: {e}")
        return [], [],  []

# 검색 버튼 클릭 시 store 업데이트 콜백 (nopedi case 전용) - legacy compatibility
@callback(
    [Output('search-selected-variety-name', 'data'),
     Output('pedigree-variety-search', 'value'),  # 검색 입력 초기화
     Output('child-nodes-store', 'data', allow_duplicate=True),  # child 노드 초기화
     Output('pedigree-path-store', 'data', allow_duplicate=True)],  # pedigree path 초기화
    Input('pedigree-search-button', 'n_clicks'),
    State('pedigree-variety-search', 'value'),
    prevent_initial_call=True
)
def handle_pedigree_search(n_clicks, selected_variety):
    """nopedi case에서 선택된 품종의 store 업데이트 (legacy compatibility)"""
    if not n_clicks or not selected_variety:
        raise PreventUpdate
    
    try:
        print(f"🔍 nopedi case pedigree search: {selected_variety}")
        # 새로운 품종 검색 시 expand 상태와 path를 초기화 (새로운 검색이므로 리셋)
        # pedigree-path-store에 검색된 품종명 추가
        path_data = [selected_variety] if selected_variety else []
        return selected_variety, None, [], path_data
        
    except Exception as e:
        print(f"❌ nopedi case pedigree search error: {e}")
        return selected_variety, None, [], []
'''


# =============================================================================
# 새로운 콜백 함수들 - Available Data 자동 스캔
# =============================================================================

# elements 변화 → Available 자동 스캔 콜백 - 중복 제거됨 (새로운 시스템에서 처리)

# 중복 callback 제거됨 - 상단 통합 callback으로 처리

# 중복 Details callback들도 제거됨 - 상단 통합 callback으로 처리

# Apply Selected(=Available) → GWAS Used 확정 콜백
'''
@callback(
    Output('gwas-state', 'data'),
    [Input('apply-selected-button', 'n_clicks')],  # Unified apply button
    [State('available-nodes-store', 'data'),
     State('gwas-state', 'data')]
)
def apply_available_to_gwas(apply_clicks, available_nodes, current_gwas_state):
    """Apply 버튼 클릭 시 available nodes 중 VCF 있는 것들을 GWAS로 로드"""
    import dash
    
    if not dash.callback_context.triggered:
        return current_gwas_state or {'status': 'idle', 'used_samples': [], 'last_error': ''}
    
    if not available_nodes:
        return {
            'status': 'error',
            'used_samples': [],
            'last_error': 'No available nodes to apply'
        }
    
    try:
        # VCF가 있는 노드들만 수집
        vcf_samples = []
        for node in available_nodes:
            if node.get('has_vcf') and node.get('vcf_id'):
                vcf_samples.append(node['vcf_id'])
        
        if not vcf_samples:
            return {
                'status': 'error',
                'used_samples': [],
                'last_error': 'No VCF data available in selected nodes'
            }
        
        print(f"🎯 GWAS 샘플 적용: {len(vcf_samples)}개 VCF 샘플")
        return {
            'status': 'ok',
            'used_samples': vcf_samples,
            'last_error': ''
        }
        
    except Exception as e:
        print(f"❌ GWAS 적용 오류: {e}")
        return {
            'status': 'error',
            'used_samples': [],
            'last_error': str(e)
        }
'''


# =============================================================================
# 🔄 새로운 Available Data 자동 스캔 시스템 - 7개 핵심 콜백
# =============================================================================

# ------------------- CALLBACK: scan_available -------------------
@callback(
    Output('available-nodes-store', 'data'),
    Input('pedigree-cytoscape', 'elements'),
    prevent_initial_call=False
)
def scan_available(elements):
    """
    elements에서 (has_vcf | has_pheno) == True 인 노드를 자동 수집.
    """
    
    return get_available_nodes_from_elements(elements)

# 중복 callback 제거됨 - 상단 통합 callback이 모든 UI를 처리
# def render_available(avail): -> 제거됨

    lines = []
    for a in avail:
        vcf = 'Available' if a.get('has_vcf') else 'Not available'
        ph  = 'Available' if a.get('has_pheno') else 'Not available'
        lines += [f"• {a['name']}", f"VCF: {vcf}", f"Phenotype: {ph}", ""]
    details = "\n".join(lines) if lines else "No data available"

    return chips, details


# ------------------- CALLBACK: apply_available -------------------
'''
@callback(
    [
        Output('applied-varieties-store', 'data'),
        Output('gwas-selected-samples-store', 'data', allow_duplicate=True),
        Output('phenotype-data-store', 'data'),
        Output('phenotype-nopedi-data-store', 'data'),   # 🆕 추가
    ],
    [
        Input('apply-selected-button', 'n_clicks'),
               # 🆕 non-pedi 데이터도 입력으로 받음
    ],
    [State('addv-apply-groups-nopedi', 'data'),
    State('available-nodes-store', 'data')],
    prevent_initial_call=True
)
def apply_available(n_clicks, group_nopedi, avail):
    """
    Pedigree + Non-Pedigree 통합 적용 처리

    - applied-varieties-store: 전체 합집합
    - gwas-selected-samples-store: has_vcf=True 인 모든 항목
    - phenotype-data-store: pedigree (P) 항목만
    - phenotype-nopedi-data-store: non-pedigree (noP) 항목만
    """
    from dash import no_update

    if not n_clicks:
        raise PreventUpdate

    nodes_p = avail or []
    nodes_np = (group_nopedi or {}).get("records", []) if isinstance(group_nopedi, dict) else []

    # 출처 태그 추가
    for n in nodes_p:
        n["source"] = "P"
    for n in nodes_np:
        n["source"] = "noP"

    all_nodes = nodes_p + nodes_np
    print(f"🔥 apply_available 호출됨: P={len(nodes_p)}, noP={len(nodes_np)}, total={len(all_nodes)}")

    # -----------------------
    # 1️⃣ GWAS 샘플
    # -----------------------
    gwas_samples = [a for a in all_nodes if a.get('has_vcf') and a.get('vcf_id')]
    vcf_ids = [str(a['vcf_id']) for a in gwas_samples]
    print(f"🧬 gwas-selected-samples-store → {len(vcf_ids)}개")

    # -----------------------
    # 2️⃣ PHENO 샘플
    # -----------------------
    def extract_it_numbers(nodes):
        its = []
        for sample in nodes:
            if not sample.get('has_pheno'):
                continue
            it_number = None
            for key in ('it_number', 'id', 'variety_id', 'IT_Number'):
                val = sample.get(key)
                if val and str(val).strip() not in ('', 'N/A', 'None'):
                    it_number = str(val).strip()
                    break
            if it_number:
                its.append(it_number)
        return list(dict.fromkeys(its))

    # 계보 기반
    it_p = extract_it_numbers(nodes_p)
    # 논계보 기반
    it_np = extract_it_numbers(nodes_np)

    pheno_store_p = [{'id': itn} for itn in it_p]
    pheno_store_np = [{'id': itn} for itn in it_np]

    print(f"📊 phenotype-data-store(P): {len(pheno_store_p)}개, phenotype-nopedi-data-store(noP): {len(pheno_store_np)}개")

    # -----------------------
    # 3️⃣ 통합 결과 반환
    # -----------------------
    return all_nodes, vcf_ids, pheno_store_p, pheno_store_np


# ------------------- CALLBACK: on_tap_expand -------------------
@callback(
    [
        Output('pedigree-cytoscape', 'elements', allow_duplicate=True),
        #Output('pedigree-cytoscape', 'stylesheet', allow_duplicate=True),
        Output('pedigree-path-store', 'data', allow_duplicate=True),
        Output('pedigree-path-child-store', 'data', allow_duplicate=True),
    ],
    Input('pedigree-cytoscape', 'tapNodeData'),
    [
        State('pedigree-cytoscape', 'elements'),
        State('pedigree-path-store', 'data'),
        State('pedigree-path-child-store', 'data'),
        State('fixed-vcf-item', 'data'),
    ],
    prevent_initial_call=True
)
def on_tap_expand(nd, base_elements, current_path, current_childs, fixed_variety):
    from dash.exceptions import PreventUpdate

    if not nd or not nd.get("id"):
        raise PreventUpdate

    nid = nd["id"]
    base_name = fixed_variety.get("processed_name") if fixed_variety else None
    current_path = current_path or []
    current_childs = current_childs or {}

    # --- base 품종 초기화
    if not any(p["type"] == "base" for p in current_path):
        if base_name:
            current_path.insert(0, {
                "index": 0, "name": base_name, "type": "base", "status": "active"
            })

    # --- 이미 있는 노드 확인
    match = next((p for p in current_path if p["name"] == nid), None)

    if match:
        # 이미 존재하는데 'removed' 상태면 복원
        if match["status"] == "removed":
            match["status"] = "active"
            print(f"♻️ 복원됨: {nid}")
        else:
            print(f"⚠️ 이미 활성 노드: {nid}")
    else:
        # 새 노드 추가
        new_index = len(current_path)
        current_path.append({
            "index": new_index, "name": nid, "type": "expand", "status": "active"
        })
        print(f"✅ 추가됨: {nid} (index={new_index})")

    # --- Cytoscape 요소 병합
    try:
        pedigree_app = get_pedigree_app()
        new_nodes, new_edges = pedigree_app.get_connected_nodes(nid, 1, 1)
        new_elements = pedigree_app.create_cytoscape_elements(new_nodes, new_edges)

        existing_ids = {e["data"]["id"] for e in base_elements if "source" not in e["data"]}
        merged = base_elements.copy()

        for ne in new_elements:
            d = ne.get("data", {})
            if "source" in d:
                if not any(
                    e.get("data", {}).get("source") == d["source"]
                    and e.get("data", {}).get("target") == d["target"]
                    for e in merged
                ):
                    merged.append(ne)
            else:
                if d.get("id") not in existing_ids:
                    merged.append(ne)

        # child-store 갱신
        current_childs[nid] = sorted(set(current_childs.get(nid, [])) | set(new_nodes) - {nid})

        available = get_available_nodes_from_elements(merged)
        
        stylesheet = build_stylesheet_with_available(
            base=[], available_ids=[x["id"] for x in available], expanded_ids=[]
        )
        
        return merged,  current_path, current_childs #stylesheet,

    except Exception as e:
        print(f"❌ expand 오류: {e}")
        raise PreventUpdate
'''




# Add dummy div for debug output (will be added to layout separately)
# html.Div(id='debug-available-nodes', style={'display': 'none'})


def get_trait_description(trait_name):
    """Trait 이름으로부터 description 등 상세 정보 가져오기"""
    try:
        if TRAIT_INFO_DF.empty:
            return {'trait_id': '', 'mapped_terms': '', 'description': '', 'source': '', 'group': ''}
        
        # Trait 정보 찾기
        trait_info = TRAIT_INFO_DF[TRAIT_INFO_DF['Trait'] == trait_name]
        
        if trait_info.empty:
            return {'trait_id': '', 'mapped_terms': '', 'description': '', 'source': '', 'group': ''}
        
        first_match = trait_info.iloc[0]
        return {
            'trait_id': first_match.get('Trait ID', ''),
            'mapped_terms': first_match.get('Mapped Terms', ''),
            'description': first_match.get('Description', ''),
            'source': first_match.get('Source', ''),
            'group': first_match.get('Group', '')
        }
        
    except Exception as e:
        print(f"❌ Trait description 조회 실패: {e}")
        return {'trait_id': '', 'mapped_terms': '', 'description': '', 'source': '', 'group': ''}


# =============================================================================
# SAMPLE SUMMARY HELPER FUNCTIONS
# =============================================================================
'''
def create_sample_summary_data(vcf_samples):
    """선택된 샘플들의 summary 데이터 생성 (Trait 단위 집계; Subtrait는 메타로만 표기)"""
    if not vcf_samples:
        return None
    try:
        gwas_result = get_gwas_data_for_samples_optimized(vcf_samples)
        if not gwas_result:
            return None

        gwas_data = gwas_result['full_data']
        if gwas_data is None or gwas_data.empty:
            return None

        sample_data_dict = gwas_result.get('sample_data', {})
        is_single = len(vcf_samples) == 1

        rows = []
        for trait in gwas_data['Trait'].dropna().unique():
            tdf = gwas_data[gwas_data['Trait'] == trait]

            subtrait = tdf['Subtrait'].iloc[0] if len(tdf) else 'Unknown'
            clean_subtrait = remove_parentheses_numbers(subtrait)

            info = get_trait_description(trait) or {}
            trait_id = info.get('trait_id', '')
            mapped_terms = info.get('mapped_terms', '')
            description = info.get('description', '')

            total_unique = int(tdf['Variant ID'].nunique())

            if is_single:
                rows.append({
                    'Trait': trait,
                    'Trait ID': trait_id,
                    'Mapped Terms': mapped_terms,
                    'Description': description,
                    'Subtrait': subtrait,
                    'Clean_Subtrait': clean_subtrait,
                    'Variant Count': total_unique
                })
            else:
                row = {
                    'Trait': trait,
                    'Trait ID': trait_id,
                    'Mapped Terms': mapped_terms,
                    'Description': description,
                    'Subtrait': subtrait,
                    'Clean_Subtrait': clean_subtrait,
                    'Total Unique Variants': total_unique
                }
                for s in vcf_samples:
                    sdf = sample_data_dict.get(s)
                    if sdf is None or sdf.empty:
                        row[f'{s}_Variant Count'] = 0
                    else:
                        c = sdf.loc[sdf['Trait'] == trait, 'Variant ID'].nunique()
                        row[f'{s}_Variant Count'] = int(c) if pd.notna(c) else 0
                rows.append(row)

        df = pd.DataFrame(rows)
        if df.empty:
            return None

        # 보기 정렬
        if is_single:
            df = df.sort_values(['Clean_Subtrait', 'Variant Count'], ascending=[True, False])
        else:
            df = df.sort_values(['Clean_Subtrait', 'Total Unique Variants'], ascending=[True, False])

        return df
    except Exception as e:
        print(f"❌ Sample summary 생성 실패: {e}")
        return None
'''


def remove_parentheses_numbers(text):
    """Subtrait에서 괄호와 숫자 제거"""
    import re
    if not text:
        return text
    # (숫자) 패턴을 찾아서 제거
    cleaned = re.sub(r'\(\d+\)', '', str(text)).strip()
    return cleaned

def get_subtrait_color(subtrait, color_map=None):
    """Subtrait에 대한 색상 반환"""
    if not color_map:
        # 기본 색상 맵 사용 (gwas_module.py의 SUBTRAIT_COLOR_MAP 참조)
        ALL_SUBTRAITS = [
            'Yield', 'Stress', 'Plantvigor', 'Biochemical', 
            'Plantgrowth_Development', 'Plant_quality', 'Biologicalprocess', 
            'Plantmorphology', 'Sterility_Fertility'
        ]
        try:
            import plotly.express as px
            ALL_SUBTRAIT_COLORS = px.colors.qualitative.Set1[:9]
            color_map = {sub: color for sub, color in zip(ALL_SUBTRAITS, ALL_SUBTRAIT_COLORS)}
            color_map['default'] = '#888888'
        except ImportError:
            # fallback colors
            color_map = {'default': '#888888'}
    
    return color_map.get(subtrait, color_map.get('default', '#888888'))

def create_sample_summary_table(summary_df, vcf_samples, primary_sample=None):
    """
    - 자동 모드: len(vcf_samples)==1 -> single columns, else multi columns
    - 버튼/토글로 들어올 때마다 최신 summary_df 기반으로 재구성
    """
    if summary_df is None or summary_df.empty:
        return html.Div([html.P("No data available", className="text-muted text-center p-3")])

    from dash import dash_table

    if not vcf_samples:
        vcf_samples = []
    if primary_sample is None:
        primary_sample = vcf_samples[0] if vcf_samples else None

    if len(vcf_samples) == 1:
        cols = ['Trait','Trait ID','Mapped Terms','Description','Subtrait','Variant Count']
        df = summary_df[[c for c in cols if c in summary_df.columns]].copy()
    else:
        base_cols = ['Trait','Trait ID','Mapped Terms','Description','Subtrait']
        dyn_cols = [c for c in summary_df.columns if c.endswith('_Variant Count')]
        cols = base_cols + dyn_cols
        df = summary_df[[c for c in cols if c in summary_df.columns]].copy()

    columns = [{'name': c, 'id': c} for c in df.columns]
    return html.Div([
        dash_table.DataTable(
            id='sample-summary-table',
            columns=columns,
            data=df.to_dict('records'),
            sort_action='native',
            filter_action='native',
            page_size=5,
            style_table={'overflowX': 'auto'},
            style_cell={'fontSize': 13, 'padding': '6px'},
        )
    ])



@callback(
    Output('multi-sample-modal', 'is_open', allow_duplicate=True),
    Input('multi-sample-toggle-btn', 'n_clicks'),
    prevent_initial_call=True
)
def open_multi_sample_modal(n_clicks):
    if n_clicks:
        return True
    raise PreventUpdate

    # 모달 닫기 콜백 (OK/Remove 버튼용)

'''
def create_all_traits_variant_presence_table(vcf_samples, min_present_count=2):
    """모든 Trait에 대해 Variant ID별 presence table 생성 + 필터 적용"""
    gwas_result = get_gwas_data_for_samples_optimized(vcf_samples)
    if not gwas_result:
        return None

    full_df = gwas_result['full_data']
    if full_df is None or full_df.empty:
        return None

    sample_data = gwas_result.get('sample_data', {})
    traits = full_df['Trait'].unique().tolist()

    rows = []
    for trait in traits:
        # trait 단위로 create_trait_variant_presence_table 로직 재사용
        trait_rows = []
        meta_df = full_df[full_df['Trait'] == trait]

        # union key 결정
        def pick_union_key(df):
            for c in ['OSAID','OSA_ID','osaid','OSA Id','OSA id','OSA ID']:
                if c in df.columns:
                    return c
            return 'Variant ID'
        union_key = pick_union_key(full_df)

        # per-sample variant set 구성
        variant_union, per_sample_sets = set(), {}
        for s in vcf_samples:
            sdf = sample_data.get(s)
            if sdf is None or sdf.empty:
                per_sample_sets[s] = set()
                continue
            xs = sdf[(sdf['Trait'] == trait)]
            ids = set(xs[union_key].astype(str).tolist())
            variant_union |= ids
            per_sample_sets[s] = ids

        # ID 매핑
        id_map = {}
        if union_key != 'Variant ID' and 'Variant ID' in meta_df.columns:
            ref = meta_df[[union_key, 'Variant ID']].drop_duplicates()
            id_map = dict(zip(ref[union_key].astype(str), ref['Variant ID'].astype(str)))

        trait_info = get_trait_description(trait)
        trait_id = trait_info.get('trait_id', '')
        mapped_terms = trait_info.get('mapped_terms', '')
        description = trait_info.get('description', '')
        subtrait = meta_df['Subtrait'].iloc[0] if len(meta_df) > 0 else 'Unknown'

        for uk in sorted(variant_union):
            display_vid = id_map.get(str(uk), str(uk))
            present_count = 0
            presence_cols = {}
            for s in vcf_samples:
                has = str(uk) in per_sample_sets.get(s, set())
                presence_cols[s] = 'Present' if has else 'Absent'
                if has: 
                    present_count += 1

            # threshold 조건 적용
            if present_count >= min_present_count:
                rows.append({
                    'Variant ID': display_vid,
                    'Trait': trait,
                    'Trait ID': trait_id,
                    'Mapped Terms': mapped_terms,
                    'Description': description,
                    'Subtrait': subtrait,
                    **presence_cols,
                    'Summary': f'{present_count}/{len(vcf_samples)}'
                })

    if not rows:
        return None

    return pd.DataFrame(rows)   

def render_snp_occurrence_table(df: pd.DataFrame, vcf_samples, table_id='snp-occurrence-table', page_size=5):
    if df is None or df.empty:
        return html.P("No GWAS data available", className="text-muted text-center p-3")

    base_cols = ['Variant ID','Trait','Trait ID','Mapped Terms','Description','Subtrait']
    sample_cols = [s for s in vcf_samples if s in df.columns]  # 존재하는 샘플 컬럼만
    extra_cols = [c for c in ['Present_Count','Total_Samples','Summary'] if c in df.columns]
    cols = [c for c in base_cols if c in df.columns] + sample_cols + extra_cols

    display_df = df[cols].copy()

    columns = [{'name': c, 'id': c} for c in display_df.columns]
    # 숫자 표현 깔끔하게 (옵션)
    for col in columns:
        if col['id'] in ['Present_Count','Total_Samples']:
            col.update({'type': 'numeric'})

    return html.Div([
        dash_table.DataTable(
            id=table_id,
            columns=columns,
            data=display_df.to_dict('records'),
            sort_action='native',
            filter_action='native',
            page_size=page_size,
            style_table={'overflowX': 'auto'},
            style_cell={'fontSize': 13, 'padding': '6px', 'textAlign': 'left'},
        )
    ])
'''

@callback(
    Output('multi-sample-modal', 'is_open'),
    Input('multi-sample-toggle-btn', 'n_clicks'),
    prevent_initial_call=True
)
def open_multi_sample_modal(n_clicks):
    if n_clicks:
        return True
    raise PreventUpdate
'''
@callback(
    Output('snp-occurrence-min-count-store', 'data'),   # ✅ 단일 소스 갱신
    [
        Input('snp-occurrence-apply-btn', 'n_clicks'),
        Input('snp-occurrence-reset-btn', 'n_clicks'),
    ],
    [
        State('snp-occurrence-min-count', 'value'),
        
        
    ],
    prevent_initial_call=True
)
def apply_or_reset_min_count(apply_clicks, reset_clicks,
                            input_value):
    if apply_clicks:
        return input_value
    elif reset_clicks:
        return 1
    return dash.no_update
'''
        

def _build_table_payload(df, vcf_samples):
    if df is None or df.empty:
        return [], []
    base_cols   = ['Variant ID','Trait','Trait ID','Mapped Terms','Description','Subtrait']
    sample_cols = [s for s in vcf_samples if s in df.columns]
    extra_cols  = [c for c in ['Present_Count','Total_Samples','Summary'] if c in df.columns]
    cols = [c for c in base_cols if c in df.columns] + sample_cols + extra_cols
    view = df[cols].copy()
    columns = [{'name': c, 'id': c} for c in view.columns]
    for col in columns:
        if col['id'] in ('Present_Count','Total_Samples'):
            col.update({'type': 'numeric'})
    return columns, view.to_dict('records')


'''
@callback(
    Output('multi-sample-modal-copy2', 'is_open'),
    Input('open-multi-sample-modal-copy2', 'n_clicks'),
    prevent_initial_call=True
)
def open_multi_sample_modal(n_clicks):
    if n_clicks:
        return True
    raise PreventUpdate
'''

def _extract_fixed_vcf(fixed_vcf_data, fixed_vcf_nopedi_data):
    """fixed-vcf-item 우선, 없으면 fixed-vcf-item-nopedi 사용. '-', 빈값 제외."""
    def pick(d):
        if isinstance(d, dict):
            val = d.get('status')
            if val is not None:
                val = str(val).strip()
                if val and val != '-':
                    return val
        return None
    return pick(fixed_vcf_data) or pick(fixed_vcf_nopedi_data)

def _dedup_preserve_order(items):
    seen, out = set(), []
    for x in items:
        if x is None:
            continue
        x = str(x).strip()
        if not x or x == '-':
            continue
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out

# ========================================
# MODULAR FILTERING SYSTEM CALLBACKS
# ========================================

# Container toggle callbacks
    
@callback(
    [Output('filter-options-collapse', 'is_open'),
    Output('filter-options-icon', 'className')],
    Input('filter-options-toggle', 'n_clicks'),
    State('filter-options-collapse', 'is_open'),
    prevent_initial_call=False
)
def toggle_filter_options(n_clicks, is_open):
    #print(n_clicks)
    """
    #v3: dbc.Collapse 방식의 Filter Options 토글
    
    Changes from #v2:
    - n_clicks 대신 is_open State 사용
    - 더 안정적인 토글 동작
    - 아이콘 자동 업데이트
    """
    if n_clicks:
        print(f"🔄 Filter Options toggle: is_open={is_open} → {not is_open}")
        new_is_open = not is_open
        icon_class = "fas fa-chevron-down" if new_is_open else "fas fa-chevron-right"
        return new_is_open, icon_class
    
    # 초기 상태
    return False, "fas fa-chevron-right"




'''
@callback(
    Output('sample-options-collapse', 'is_open'),
    Output('sample-options-icon', 'className'),
    Input('sample-options-toggle', 'n_clicks'),
    State('sample-options-collapse', 'is_open'),
    prevent_initial_call=True
)
def toggle_sample_options(n, is_open):
    if n:
        new_state = not is_open
        icon = "fas fa-chevron-down" if new_state else "fas fa-chevron-right"
        return new_state, icon
    return is_open, "fas fa-chevron-right"

@callback(
    Output('individual-sample-section', 'style'),
    Output('merge-sample-section', 'style'),
    Output('unique-sample-section', 'style'),
    Input('sample-mode-radio', 'value')
)
def update_sample_mode(selected_mode):
    print("🔄 selected_mode:", selected_mode)

    base = {
        'marginTop': 10, 
        'padding': '10px', 
        'border': '1px solid #ccc'
    }
    hidden = {**base, 'display': 'none'}
    visible = {**base, 'display': 'block'}

    if selected_mode == 'each':
        return visible, hidden, hidden
    elif selected_mode == 'merge':
        return hidden, visible, hidden
    elif selected_mode == 'unique':
        return hidden, hidden, visible
    return no_update, no_update, no_update
'''
@callback(
    Output("unique-mode-enabled", "value",allow_duplicate=True),          # 스위치 자동 off 제어
    Output("unique-mode-warning", "children",allow_duplicate=True),       # 경고 메시지
    Output("sample-mode-store", "data",allow_duplicate=True),             # 모드 업데이트
    Input("gwas-selected-samples-store", "data"),    # 샘플 개수 변화 감시
    State("unique-mode-enabled", "value"),           # 현재 unique 스위치 상태
    prevent_initial_call=True
)
def update_sample_mode(selected_samples, unique_enabled):
    """
    - unique_enabled: True → 스위치 ON
                      False → 스위치 OFF
                      None → 스위치가 아직 layout에 없음
    """
    n = len(selected_samples) if selected_samples else 0

    # 스위치가 아직 생기지 않은 상태면 — 아무 작업 안 함
    if unique_enabled is None:
        return no_update, no_update, no_update

    # 스위치가 꺼져 있을 때
    if not unique_enabled:
        if n == 1:
            return False, "", "each"
        elif n >= 2:
            return False, "", "merge"
        else:
            return False, "", no_update

    # 스위치가 켜져 있을 때
    if 2 <= n <= 5:
        return True, "", "unique"

    warning = f"❌ Unique mode requires 2–5 samples (currently {n})"
    if n == 1:
        return False, warning, "each"
    elif n >= 6:
        return False, warning, "merge"

    return False, warning, no_update



@callback(
    [
        Output('trait-section', 'style'),
        Output('pvalue-section', 'style'),
        Output('maf-section', 'style'),
        Output('snp-presence-section', 'style'),
        Output('unique-section', 'style')  # ✅ 추가됨
    ],
    Input('filter-radio', 'value')
)
def toggle_filter_sections(selected_filters):
    base_style = {
        'padding': '10px',
        'border': '1px solid #ccc',
        'marginBottom': '10px',
        'borderRadius': '4px',
        'backgroundColor': '#f8f9fa'
    }
    hidden_style = {**base_style, 'display': 'none'}
    visible_style = {**base_style, 'display': 'block'}

    selected_filters = selected_filters or []

    trait_style = visible_style if 'trait' in selected_filters else hidden_style
    pvalue_style = visible_style if 'pvalue' in selected_filters else hidden_style
    maf_style = visible_style if 'maf' in selected_filters else hidden_style
    snp_presence_style = visible_style if 'snp_presence' in selected_filters else hidden_style
    unique_style = visible_style if 'unique' in selected_filters else hidden_style  # ✅ 추가됨

    return trait_style, pvalue_style, maf_style, snp_presence_style, unique_style

@callback(
    Output("unique-mode-enabled", "value",allow_duplicate=True),          # 스위치 on/off 제어
    Output("unique-mode-warning", "children",allow_duplicate=True),       # 경고 메시지
    Output("sample-mode-store", "data",allow_duplicate=True),             # 현재 모드(each / merge / unique)
    Input("unique-mode-enabled", "value"),           # 스위치 클릭 감지
    State("gwas-selected-samples-store", "data"),    # 선택된 샘플 목록
    prevent_initial_call=True
)
def update_unique_filter_state(enabled, selected_samples):
    """
    Unique mode 스위치 상태, 경고 메시지, 샘플 모드(each/merge/unique)를 동시에 관리.
    - 2~5개의 샘플에서만 unique 모드 허용
    - 조건 벗어나면 스위치를 자동 off하고 경고 표시
    - 스위치 off 시에도 샘플 수에 따라 each/merge 모드 유지
    """
    n = len(selected_samples) if selected_samples else 0

    # --- 1️⃣ Unique 스위치가 ON인 경우 ---
    if enabled:
        if 2 <= n <= 5:
            # 정상 범위 → unique 모드 유지
            return True, "", "unique"

        # 범위 벗어남 → 스위치 off + 경고
        warning = f"❌ Unique mode requires 2–5 samples (currently {n})"
        if n == 1:
            return False, warning, "each"
        elif n >= 6:
            return False, warning, "merge"
        else:
            # n == 0 인 경우 등 → 기존 모드 유지
            return False, warning, no_update

    # --- 2️⃣ Unique 스위치가 OFF인 경우 ---
    else:
        if not selected_samples:
            # 샘플 없으면 아무것도 변경하지 않음
            return False, "", no_update

        if n == 1:
            return False, "", "each"
        elif n >= 2:
            return False, "", "merge"

        return False, "", no_update







# Connect subtrait selection to the plot container (using existing integrated-trait-bar-container logic)

        
@callback(
    Output('integrated-selected-traits-store', 'data', allow_duplicate=True),
    Input('integrated-trait-barplot', 'clickData'),
    State('integrated-selected-traits-store', 'data'),
    prevent_initial_call=True
)
def handle_integrated_trait_barplot_click(clickData, current_traits):
    """
    Handle barplot click → toggle selected traits
    """
    current_traits = current_traits or []

    if not clickData or not clickData.get('points'):
        return current_traits

    # ✅ x 값이 아닌 customdata에서 Trait 이름 가져오기
    clicked_trait = clickData['points'][0].get('customdata')

    if not clicked_trait:
        print("⚠️ No trait name found in clickData")
        return current_traits

    # Toggle
    if clicked_trait in current_traits:
        current_traits.remove(clicked_trait)
        print(f"🔄 Trait '{clicked_trait}' removed from selection")
    else:
        current_traits.append(clicked_trait)
        print(f"✅ Trait '{clicked_trait}' added to selection")

    return current_traits
  



    

    '''
    @callback(
        Output('integrated-trait-barplot', 'figure'),
        [
            Input('integrated-selected-traits-store', 'data')
        ],
        State('copy2gwas-selected-samples-store', 'data'),
        State('integrated-subtrait-dropdown', 'value'),
        prevent_initial_call=True
    )
    def update_barplot_on_trait_selection(selected_traits, selected_samples, selected_subtraits):
        if not selected_samples:
            raise PreventUpdate
        
        summary_df = create_sample_summary_data(selected_samples)
        if summary_df is None or summary_df.empty:
            raise PreventUpdate
        
        # 필터링
        if selected_subtraits:
            df_plot = summary_df[summary_df['Subtrait'].isin(selected_subtraits)].copy()
        else:
            df_plot = summary_df.copy()
        
        return create_enhanced_clickable_barplot(df_plot, selected_samples, selected_traits)
    '''
def create_sample_summary_data_from_combo(
    combo_json, vcf_samples,
    maf_cutoff=None, maf_enabled=False,
    presence_threshold=None, presence_enabled=False
):
    if not combo_json or combo_json.get("status") != "ready":
        return None, None, {}

    filter_meta = {"maf": {}, "presence": {}}

    # --- Base variant_only ---
    df_vo = _records_to_df(combo_json.get("variant_only"))
    df_vp = _records_to_df(combo_json.get("variant_pvalue"))  # ✅ 항상 로드
    if df_vo.empty:
        return None, df_vp, {}
    if not df_vp.empty and "Subtrait" in df_vp.columns:
        df_vp["Clean_Subtrait"] = df_vp["Subtrait"].apply(remove_parentheses_numbers)

    # --- Case 1: MAF filter ---
    if maf_enabled and maf_cutoff is not None and not df_vp.empty:
        df_vp = df_vp[df_vp["MAF"].astype(float) >= maf_cutoff]
        valid_variants = set(df_vp["Variant ID"].unique())
        df_vo = df_vo[df_vo["Variant ID"].isin(valid_variants)]

        filter_meta["maf"] = {
            "enabled": True,
            "cutoff": maf_cutoff,
            "variant_ids": list(valid_variants)
        }
    else:
        filter_meta["maf"] = {"enabled": False}

    # --- Case 2: Presence filter ---
    if presence_enabled and presence_threshold is not None and vcf_samples:
        gt_cols = [f"{s}_GT" for s in vcf_samples if f"{s}_GT" in df_vo.columns]
        if gt_cols:
            df_vo["_present_count"] = df_vo[gt_cols].notna().sum(axis=1)
            valid_df = df_vo[df_vo["_present_count"] >= presence_threshold]
            valid_variants = set(valid_df["Variant ID"].unique())
            df_vo = valid_df.drop(columns=["_present_count"])

            # df_vp 도 동기화 (있을 경우)
            if not df_vp.empty:
                df_vp = df_vp[df_vp["Variant ID"].isin(valid_variants)]

            filter_meta["presence"] = {
                "enabled": True,
                "threshold": presence_threshold,
                "total_samples": len(gt_cols),
                "variant_ids": list(valid_variants)
            }
        else:
            filter_meta["presence"] = {"enabled": False}
    else:
        filter_meta["presence"] = {"enabled": False}

    # --- Summary DF 생성 ---
    if df_vo.empty:
        return None, df_vp, filter_meta

    rows = []
    is_single = len(vcf_samples) == 1
    for trait in df_vo["Trait"].dropna().unique():
        tdf = df_vo[df_vo["Trait"] == trait]
        subtrait = tdf["Subtrait"].iloc[0] if "Subtrait" in tdf.columns and len(tdf) else "Unknown"
        clean_subtrait = remove_parentheses_numbers(subtrait)

        info = get_trait_description(trait) or {}
        trait_id = info.get("trait_id", "")
        mapped_terms = info.get("mapped_terms", "")
        description = info.get("description", "")

        total_unique = int(tdf["Variant ID"].nunique())

        if is_single:
            rows.append({
                "Trait": trait,
                "Trait ID": trait_id,
                "Mapped Terms": mapped_terms,
                "Description": description,
                "Subtrait": subtrait,
                "Clean_Subtrait": clean_subtrait,
                "Variant Count": total_unique
            })
        else:
            row = {
                "Trait": trait,
                "Trait ID": trait_id,
                "Mapped Terms": mapped_terms,
                "Description": description,
                "Subtrait": subtrait,
                "Clean_Subtrait": clean_subtrait,
                "Total Unique Variants": total_unique
            }
            for s in vcf_samples:
                gt_col = f"{s}_GT"
                if gt_col in tdf.columns:
                    cnt = (tdf[gt_col].notna() & (tdf[gt_col].astype(str) != "")).sum()
                    row[f"{s}_Variant Count"] = int(cnt)
                else:
                    row[f"{s}_Variant Count"] = 0
            rows.append(row)

    df = pd.DataFrame(rows)
    if df.empty:
        return None, df_vp, filter_meta

    # 정렬
    if is_single:
        df = df.sort_values(["Clean_Subtrait", "Variant Count"], ascending=[True, False])
    else:
        df = df.sort_values(["Clean_Subtrait", "Total Unique Variants"], ascending=[True, False])

    return df, df_vp, filter_meta

def build_subtrait_options_from_combo(combo_json):
    """
    combo_json['variant_pvalue']에서 Subtrait 옵션 생성
    - df['Subtrait'].unique() 기반
    - label: 괄호 제거된 순수 이름
    - value: 원래 Subtrait 문자열 그대로 (예: 'Plantvigor(23)')
    """
    if not combo_json or combo_json.get('status') != 'ready':
        return []

    df = _records_to_df(combo_json.get('variant_pvalue'))
    if df.empty or 'Subtrait' not in df.columns:
        return []

    # Subtrait 리스트 뽑기
    st_list = df['Subtrait'].dropna().astype(str).unique().tolist()
    st_list.sort()

    options = []
    for st in st_list:
        # label에서는 괄호 포함된 숫자 부분 제거
        clean_label = st.split("(")[0].strip().replace('_', ' ')
        options.append({
            "label": clean_label.strip(),
            "value": st  # 원래 값 그대로 유지
        })

    return options

@callback(
    Output("download-variant-only", "data"),
    Input("btn-download-variant-only", "n_clicks"),
    State("gwas-combo-store", "data"),
    prevent_initial_call=True
)
def download_variant_only(n_clicks, combo_json):
    if not combo_json or combo_json.get("status") != "ready":
        return None

    df_vo = _records_to_df(combo_json.get("variant_only"))
    if df_vo.empty:
        return None

    # CSV 파일로 반환
    return dcc.send_data_frame(df_vo.to_csv, "variant_only.csv", index=False)



@callback(
    Output("maf-filter-store", "data"),
    [Input("maf-enabled", "value"),
     Input("maf-cutoff", "value")],
    prevent_initial_call=True
)
def update_maf_filter_state(enabled, cutoff):
    """MAF filter 상태를 Store에 저장"""
    state = {
        "enabled": bool(enabled),
        "cutoff": float(cutoff) if cutoff is not None else 0.05
    }
    print(f"🎯 MAF filter state updated: {state}")
    return state   


@callback(
    Output('snp-occurrence-store', 'data'),
    Output('snp-presence-status', 'children'),
    [Input('snp-occurrence-enabled', 'value'),
     Input('snp-occurrence-min-count', 'value')],
    prevent_initial_call=True
)
def update_presence_store(enabled, min_count):
    if not enabled:
        return {"enabled": False, "threshold": None}, "❌ SNP Presence filter disabled"

    threshold = int(min_count) if min_count else 1
    status_text = f"✅ Active: Variants present in ≥ {threshold} samples"

    return {"enabled": True, "threshold": threshold}, status_text


def create_selected_traits_table(selected_traits, trait_table):
    """선택된 traits에 대한 상세 정보 테이블 생성"""
    import dash_table

    if not selected_traits or trait_table.empty:
        return html.Div()

    filtered_table = trait_table[trait_table['Trait'].isin(selected_traits)].copy()
    if filtered_table.empty:
        return html.P("No matching trait info found.", style={'color': '#666'})

    # Rename Group → Subtrait
    filtered_table = filtered_table.rename(columns={'Group': 'Subtrait'})

    display_columns = ['Trait', 'Trait ID', 'Mapped Terms', 'Description', 'Subtrait']
    table_data = filtered_table[display_columns].to_dict('records')

    return dash_table.DataTable(
        data=table_data,
        columns=[{"name": col, "id": col} for col in display_columns],
        style_cell={
            'textAlign': 'left',
            'padding': '8px',
            'fontSize': '12px',
            'fontFamily': 'Arial, sans-serif',
            'overflow': 'hidden',
            'textOverflow': 'ellipsis'
        },
        style_header={
            'backgroundColor': 'rgb(240,240,240)',
            'fontWeight': 'bold',
            'textAlign': 'center'
        },
        style_data={'whiteSpace': 'normal', 'height': 'auto'},
        style_table={'overflowX': 'auto'},
        tooltip_data=[
            {
                column: {'value': str(row[column]), 'type': 'markdown'}
                for column in display_columns
            }
            for row in table_data
        ],
        tooltip_duration=None,
        page_size=10,
        sort_action="native"
    )



@callback(
    Output("final-df-store", "data"),
    Input("sample-mode-store", "data"),               # ✅ 변경됨
    Input("gwas-selected-samples-store", "data"),     # ✅ 새로 추가됨
    Input("processed-df-store", "data"),              # ✅ 기존 유지
    prevent_initial_call=True
)
def update_final_df_and_size(mode, selected_samples, processed):
    """
    최종 분석용 데이터(final-df-store)를 구성합니다.
    - mode: "each", "merge", "unique"
    - selected_samples: 현재 선택된 GWAS 샘플 리스트
    - processed: preprocess_df_vp()에서 만든 처리 결과
    """
    import pandas as pd

    if not processed or not mode or not selected_samples:
        return no_update

    # ─────────────────────────────
    # EACH 모드 (샘플 1개)
    # ─────────────────────────────
    if mode == "each":  
        sample = selected_samples[0]
        if sample in processed["each"]:
            df = pd.DataFrame(processed["each"][sample])
            print(f"✅ EACH mode: {sample} ({len(df)} rows)")
            return processed["each"][sample]
        print(f"⚠️ EACH mode sample not found in processed data: {sample}")
        return None

    # ─────────────────────────────
    # MERGE 모드 (샘플 ≥2개)
    # ─────────────────────────────
    elif mode == "merge":
        df = pd.DataFrame(processed["merge"])
        print(f"✅ MERGE mode: {len(selected_samples)} samples merged ({len(df)} rows)")
        return processed["merge"]

    # ─────────────────────────────
    # UNIQUE 모드
    # ─────────────────────────────
    elif mode == "unique":
        if isinstance(processed["unique"], dict):
            selected_unique = {
                s: processed["unique"].get(s, [])
                for s in selected_samples if s in processed["unique"]
            }
            total_rows = sum(len(v) for v in selected_unique.values())
            #print(f"✅ UNIQUE mode: {len(selected_samples)} samples, total {len(total_rows)} unique rows")
            return selected_unique
        else:
            print(f"⚠️ UNIQUE data is not dict type: {type(processed['unique'])}")
            return processed["unique"]

    print(f"⚠️ Unknown mode: {mode}")
    return None

#Output('gwas-plot-mode-radio-container', 'style') 제거
@callback(
    [Output('gwas-sample-scatter', 'figure', allow_duplicate=True),
    Output('gwas-sample-scatter-group-data', 'data', allow_duplicate=True),
     Output('gwas-sample-scatter-container', 'style',allow_duplicate=True),
     ],
    [Input('final-df-store', 'data'),
     Input('sample-mode-store', 'data'),
     Input('pvalue-cutoff', 'value'),
     Input('integrated-selected-traits-store', 'data'),
      Input('gwas-click-store', 'data'),   
     #Input('scatter-refresh-trigger', 'data'),
     
     ],
    prevent_initial_call=True
)
def update_integrated_scatter(final_df, mode, pvalue_cutoff,selected_traits, click_store):#,selected_traits
    import pandas as pd, numpy as np
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    from dash import ctx
    trig = ctx.triggered_id


    if not final_df or not mode:
        return go.Figure().add_annotation(
            text="No data available",
            x=0.5, y=0.5, xref="paper", yref="paper",
            showarrow=False, font=dict(size=16, color="gray")
        ), no_update, no_update
    
        
    

    if trig == 'gwas-click-store':
        # None -> 패스, []만 재생성, 길이>0 -> 패스
        if click_store is None or len(click_store) != 0:
            raise PreventUpdate
    


    # ---------------- 데이터 준비 ----------------
    if mode == "each":
        df = pd.DataFrame(final_df)
        gt_cols = [c for c in df.columns if c.endswith("_GT")]
        if gt_cols:
            df["GT"] = df[gt_cols[0]]

    elif mode == "merge":
        df = pd.DataFrame(final_df)
        gt_cols = [c for c in df.columns if c.endswith("_GT")]
        df["GTs"] = df[gt_cols].apply(
            lambda row: {
                col.replace("_GT", ""): row[col]
                for col in gt_cols if pd.notna(row[col])
            },
            axis=1
        ).astype(str)

    elif mode == "unique":
        records = []
        for sample, rows in final_df.items():
            for r in rows:
                r["Sample"] = sample
                records.append(r)
        df = pd.DataFrame(records)
        gt_cols = [f"{s}_GT" for s in final_df.keys()]

    else:
        return go.Figure().add_annotation(
            text="Unknown mode",
            x=0.5, y=0.5, xref="paper", yref="paper",
            showarrow=False, font=dict(size=16, color="red")
        ), no_update, no_update

    if df.empty:
        return go.Figure().add_annotation(
            text="Empty dataset",
            x=0.5, y=0.5, xref="paper", yref="paper",
            showarrow=False, font=dict(size=16, color="gray")
        ), no_update, no_update

    # ---------------- P-value 처리 ----------------
    df["P-value"] = pd.to_numeric(df["P-value"], errors="coerce")
    df["neg_log10_pvalue"] = df["P-value"].apply(lambda x: -np.log10(x) if pd.notna(x) and x > 0 else None)

    # ---------------- Subtrait 색상 매핑 ----------------
    df["Color"] = df["Subtrait"].map(lambda s: SUBTRAIT_COLOR_MAP.get(s, "#888888"))
    #df["Clean_Subtrait"] = df["Subtrait"].astype(str).str.replace(r"\(\d+\)", "", regex=True).str.strip()

    # ---------------- Subplot 구조 ----------------
    chromosomes = sorted(df["Chromosome"].dropna().unique(), key=lambda x: int(x))
    cols = min(4, len(chromosomes))
    rows = (len(chromosomes) + cols - 1) // cols
    subplot_titles = [f"CHR{c}" for c in chromosomes]

    fig = make_subplots(
        rows=rows, cols=cols,
        subplot_titles=subplot_titles,
        vertical_spacing=0.12,
        horizontal_spacing=0.08
    )

    marker_symbols = ["circle", "square", "diamond", "cross", "x", "triangle-up", "star"]

    # 전역 legend tracker (중복 방지)
    seen_legends = set()

    # ---------------- Chromosome 루프 ----------------
    for i, chr_id in enumerate(chromosomes):
        row, col = i // cols + 1, i % cols + 1
        subdf = df[df["Chromosome"] == chr_id]
        if subdf.empty:
            continue

        if mode in ["each", "merge"]:
            # Subtrait 단위
            for subtrait, grp in subdf.groupby("Clean_Subtrait"):
                legend_name = subtrait
                legend_group = f"subtrait_{subtrait}"
                showlegend = legend_name not in seen_legends
                line_width = 1.3 if grp["Trait"].iloc[0] in (selected_traits or []) else 0.5
                gt_column = grp["GTs"] if mode == "merge" else grp["GT"]
                gt_label = "GTs" if mode == "merge" else "GT"
                fig.add_trace(go.Scattergl(
                    x=grp["Position"], y=grp["neg_log10_pvalue"],
                    mode="markers",
                    marker=dict(color=grp["Color"].iloc[0], size=7,
                                line=dict(width=line_width , color="black")),
                    name=legend_name,
                    legendgroup=legend_group,
                    showlegend=showlegend,
                    customdata=np.stack([
                        grp["Chromosome"],
                        grp["P-value"],
                        grp["Clean_Subtrait"],
                        grp["Minor Allele"],
                        grp['Variant ID'],
                        gt_column
                    ], axis=-1),
                    hovertemplate=(
                        "Chr: %{customdata[0]}<br>"
                        "Pos: %{x}<br>"
                        "P: %{customdata[1]:.2e}<br>"
                        "-log₁₀(P): %{y:.2f}<br>"
                        "Subtrait: %{customdata[2]}<br>"
                        "Minor Allele: %{customdata[3]}<br>"
                        "Variant ID: %{customdata[4]}<br>"
                         f"{gt_label}: "+"%{customdata[5]}<extra></extra>"
                    )
                ), row=row, col=col)

                seen_legends.add(legend_name)

        elif mode == "unique":
            n_samples = len(final_df.keys()) if isinstance(final_df, dict) else 0
            total_records = sum(len(rows) for rows in final_df.values())
            if n_samples > 4:
                return go.Figure().add_annotation(
                    text="❌ UNIQUE mode supports ≤4 samples only.\nScatter plot is not available.",
                    x=0.5, y=0.5, xref="paper", yref="paper",
                    showarrow=False, font=dict(size=16, color="red")
                ), no_update, no_update
            elif total_records == 0:
                return go.Figure().add_annotation(
                    text="⚠ No unique variants available.",
                    x=0.5, y=0.5, xref="paper", yref="paper",
                    showarrow=False, font=dict(size=16, color="gray")
                ), no_update, no_update
            else:
                
                for j, sample in enumerate(subdf["Sample"].unique()):
                    sample_df = subdf[subdf["Sample"] == sample]
                    # Subtrait 정렬 순서
                    for subtrait in sorted(sample_df["Clean_Subtrait"].dropna().unique()):
                        grp = sample_df[sample_df["Clean_Subtrait"] == subtrait]
                        legend_name = f"{sample} | {subtrait}"
                        legend_group = f"{sample}_{subtrait}"
                        showlegend = legend_name not in seen_legends
                        line_width = 1.3 if grp["Trait"].iloc[0] in (selected_traits or []) else 0.5

                        fig.add_trace(go.Scattergl(
                            x=grp["Position"], y=grp["neg_log10_pvalue"],
                            mode="markers",
                            marker=dict(
                                color=grp["Color"].iloc[0],
                                symbol=marker_symbols[j % len(marker_symbols)],
                                size=7, line=dict(width=line_width, color="black")
                            ),
                            name=legend_name,
                            legendgroup=legend_group,
                            showlegend=showlegend,
                            customdata=np.stack([
                                grp["Chromosome"],
                                grp["P-value"],
                                grp["Clean_Subtrait"],
                                grp["Minor Allele"],
                                grp['Variant ID'],
                                grp.get(f"{sample}_GT", [""]*len(grp))
                            ], axis=-1),
                            hovertemplate=(
                                f"<b>Sample:</b> {sample}<br>"
                                "Chr: %{customdata[0]}<br>"
                                "Pos: %{x}<br>"
                                "P: %{customdata[1]:.2e}<br>"
                                "-log₁₀(P): %{y:.2f}<br>"
                                "Subtrait: %{customdata[2]}<br>"
                                "Minor Allele: %{customdata[3]}<br>"
                                "Variant ID: %{customdata[4]}<br>"
                                "GT: %{customdata[5]}<extra></extra>"
                            )
                        ), row=row, col=col)

                        seen_legends.add(legend_name)

        # cutoff 선
        if pvalue_cutoff and pvalue_cutoff > 0:
            fig.add_hline(
                y=pvalue_cutoff, line_dash="dash",
                line_color="red", line_width=2,
                row=row, col=col
            )

    # ---------------- 레이아웃 ----------------
    fig.update_layout(
        height=800,
        showlegend=True,
        margin=dict(t=50, b=40, l=50, r=20),
        xaxis_title="Position",
        yaxis_title="-log₁₀(P)"
    )
    used_samples = [c.replace("_GT", "") for c in gt_cols]
    return fig, {"used_samples": used_samples}, {'display': 'block'}
'''
@callback(
    [Output('gwas-sample-scatter-wrapper', 'style',allow_duplicate=True),
     Output('gwas-sample-table-wrapper', 'style',allow_duplicate=True)],
    Input('gwas-plot-mode-radio', 'value'),
    prevent_initial_call=True
)
def toggle_scatter_table(tab_value):
    if tab_value == "scatter":
        return {'display': 'block'}, {'display': 'none'}
    elif tab_value == "table":
        return {'display': 'none'}, {'display': 'block'}
    else:
        return {'display': 'none'}, {'display': 'none'}
        '''

@callback(
    Output('group-options-collapse', 'is_open'),
    Output('group-options-icon', 'className'),
    Input('group-options-toggle', 'n_clicks'),
    State('group-options-collapse', 'is_open'),
    prevent_initial_call=True
)
def toggle_group_section(n_clicks, is_open):
    if n_clicks:
        new_state = not is_open
        icon = "fas fa-chevron-down" if new_state else "fas fa-chevron-right"
        return new_state, icon
    return is_open, no_update

@callback(
    Output('trait-group-container', 'style'),
    Output('variant-group-container', 'style'),
    Output('sample-group-container', 'style'),
    Output('copy2selected-traits-store', 'data',allow_duplicate=True),
    Input('group-type-radio', 'value'),
    State('integrated-selected-traits-store', 'data'),
    prevent_initial_call=True
)
def switch_group_container(selected_type, integrated_traits):
    """
    1. group-type-radio 값에 따라 group container 표시 전환
    2. group-type이 'trait'일 때만 integrated-selected-traits-store → copy2selected-traits-store 복사
    """
    hidden = {'display': 'none'}
    visible = {'display': 'block'}
    integrated_traits = integrated_traits or []

    if selected_type == 'trait':
        # ✅ container 표시 + 통합 store → copy2selected-traits-store 복사
        return visible, hidden, hidden, integrated_traits

    elif selected_type == 'variant':
        # variant 활성화 시에는 복사하지 않음
        return hidden, visible, hidden, dash.no_update

    elif selected_type == 'sample':
        # sample 활성화 시에도 복사하지 않음
        return hidden, hidden, visible, dash.no_update

    # 기본값 (fallback)
    return visible, hidden, hidden, dash.no_update

@callback(
    [Output({'type': 'group-block', 'group_type': MATCH, 'index': i}, 'style') for i in range(3, 6)],
    Input({'type': 'add-group-btn', 'group_type': MATCH}, 'n_clicks'),
    [Input({'type': 'remove-group-btn', 'group_type': MATCH, 'index': i}, 'n_clicks') for i in range(3, 6)],
    [State({'type': 'group-block', 'group_type': MATCH, 'index': i}, 'style') for i in range(3, 6)],
    prevent_initial_call=True
)
def update_groups(add_click, r3, r4, r5, g3, g4, g5):
    states = [g3, g4, g5]
    trigger = ctx.triggered_id

    if not trigger:
        return states

    # ✅ Add
    if trigger["type"] == "add-group-btn":
        for idx, s in enumerate(states):
            if s["display"] == "none":
                new_state = s.copy()
                new_state["display"] = "flex"   # block 대신 flex
                states[idx] = new_state
                return states
        return states

    # ✅ Remove
    if trigger["type"] == "remove-group-btn":
        remove_idx = trigger["index"] - 3
        states[remove_idx] = {**states[remove_idx], "display": "none"}
        return states

    return states

@callback(
    Output("trait-group-container", "children"),
    Input("integrated-selected-traits-store", "data")
)
def render_trait_grouping(selected_traits):
    if not selected_traits or len(selected_traits) < 3:
        return html.Div(
            "⚠️ Please select at least 3 traits to enable grouping.",
            style={'color': 'red', 'fontSize': '13px', 'padding': '6px'}
        )

    return html.Div([
        # --- 상단 설명 header ---
        html.Div([
            html.I(className="fas fa-info-circle", style={'color': '#007bff', 'marginRight': '6px'}),
            html.Span("Trait grouping", style={
                'fontWeight': 'bold',
                'color': '#007bff',
                'marginRight': '6px'
            }),
            html.Span(
                "Available only when ≥ 3 traits are selected.",
                style={'fontSize': '12px', 'color': '#555'}
            )
        ], style={
            'display': 'flex',
            'alignItems': 'center',
            'marginBottom': '6px'
        }),

        # --- Selected Traits Section (source summary) ---
        html.Div(
            id="selected-traits-section",
            children=[
                html.Div(
                    [
                        html.Span("Selected Traits:", 
                                  style={'fontWeight': 'bold', "color": "#007bff",
                                         "marginRight": "6px", 'fontSize': '14px'}),
                        dbc.Button(
                            [html.I(id="selected-traits-icon", className="fas fa-table"),
                             html.Span("(closed)", id="selected-traits-status",
                                       style={'marginLeft': '5px', 'fontSize': '12px'})],
                            id="selected-traits-toggle",
                            color="link",
                            size="sm",
                            style={
                                'marginLeft': '8px',
                                'cursor': 'pointer',
                                'fontSize': '14px',
                                'textDecoration': 'none',
                                'border': 'none',
                                'padding': 0
                            }
                        )
                    ],
                    style={'display': 'flex', 'alignItems': 'center'}
                ),

                html.Pre(
                    id="selected-traits-summary",
                    style={'margin': 0, 'whiteSpace': 'pre-wrap', 'fontSize': '12px',
                           'textAlign': 'left'}
                ),

                html.Div(
                    id="integrated-selected-traits-table",
                        style={
                        'marginTop': '10px',
                        'maxHeight': '0px',
                        'overflow': 'hidden',
                        'transition': 'max-height 0.3s ease-out'
                    }
                ),
            ],
            style={
                'border': '1px solid #ccc',
                'borderRadius': '4px',
                'padding': '8px',
                'width': '100%',
                'marginTop': '6px',
                'backgroundColor': '#f8f9fa'
            }
        ),

        # --- Add Group 버튼 ---
        dbc.Button(
            "➕ Add Group",
            id={'type': 'add-group-btn', 'group_type': 'trait'},
            color="secondary",
            size="sm",
            style={
                'padding': '2px 6px',
                'fontSize': '12px',
                'marginTop': '10px',
                'marginBottom': '8px'
            }
        ),

        # --- 구분선 ---
        html.Hr(style={
            'marginTop': '10px',
            'marginBottom': '10px',
            'borderTop': '1px solid #ddd'
        }),

        # --- Group block 리스트 ---
        html.Div([
            create_group_block(1, visible=True, removable=False, group_type="trait"),
            create_group_block(2, visible=True, removable=False, group_type="trait"),
            create_group_block(3, visible=False, group_type="trait"),
            create_group_block(4, visible=False, group_type="trait"),
            create_group_block(5, visible=False, group_type="trait"),
        ], style={
            'display': 'flex',
            'flexWrap': 'wrap',
            'gap': '10px',
            'justifyContent': 'flex-start',
            'padding': '5px'
        }),
        html.Hr(style={'marginTop': '3px', 'marginBottom': '3px'}),
            html.Div(
                id="group-enable-container",
                style={
                    "display": "flex",
                    "flexDirection": "column",
                    "alignItems": "flex-end",
                    "margin": "10px 15px 5px 10px",
                    "gap": "4px"
                },
                children=[
                    dbc.Button(
                        "Apply Grouping to Plots",
                        id="group-enable-switch",
                        n_clicks=0,
                        color="secondary",
                        size="sm",
                        style={"width": "220px"}
                    )
                ]
            )

    ],
    style={
        'border': '1px solid #ccc',
        'borderRadius': '6px',
        'padding': '12px 15px',
        'backgroundColor': '#ffffff',
        'marginBottom': '10px'
    })


@callback(
    Output({'type': 'group-dropdown', 'group_type': 'trait', 'index': ALL}, 'options'),
    Output({'type': 'group-dropdown', 'group_type': 'trait', 'index': ALL}, 'value'),
    Output("presence-sets-store", "data", allow_duplicate=True),
    Input("integrated-selected-traits-store", "data"),
    Input({'type': 'group-dropdown', 'group_type': 'trait', 'index': ALL}, 'value'),
    prevent_initial_call=True
)
def update_trait_group_options(selected_traits, group_values):
    if not selected_traits or len(selected_traits) < 3:
        return [[] for _ in group_values], [None for _ in group_values], None

    consumed = set([t for vals in group_values if vals for t in vals])
    options_list, values_out = [], []

    for vals in group_values:
        available = [t for t in selected_traits if (t not in consumed or t in (vals or []))]
        opts = [{'label': t, 'value': t} for t in available]
        options_list.append(opts)

        # ✅ 기존 선택값 중 유효한 것만 유지
        if vals:
            valid_vals = [v for v in vals if v in available]
            values_out.append(valid_vals if valid_vals else None)
        else:
            values_out.append(None)

    return options_list, values_out, None

@callback(
    Output("variant-group-container", "children"),
    Input("filter-meta-store", "data")
)
def render_variant_grouping(filter_meta):
    #print(filter_meta)
    #print(len(filter_meta))
    presence_sets = filter_meta.get("presence_sets", []) if filter_meta else []
    #print(presence_sets)
    if not presence_sets or len(presence_sets) < 2:
        return html.Div(
            "⚠️ Please define at least 2 distinct SNP Presence Threshold sets to enable variant grouping.",
            style={'color': 'red', 'fontSize': '13px', 'padding': '6px'}
        )

    return html.Div([
        # --- 상단 설명 + Add 버튼 ---
        html.Div([
            html.Div([
                html.I(className="fas fa-info-circle", style={'color': '#007bff', 'marginRight': '6px'}),
                html.Span("Variant ID grouping", style={
                    'fontWeight': 'bold',
                    'color': '#007bff',
                    'marginRight': '6px'
                }),
                html.Span(
                    "Requires ≥2 distinct SNP Presence Threshold sets.",
                    style={'fontSize': '12px', 'color': '#555'}
                )
            ], style={'marginBottom': '6px', 'display': 'flex', 'alignItems': 'center'}),

            html.Div([
                html.Ul([
                    html.Li("Set thr=5 → variants present in ≥5 samples (unique group)"),
                    html.Li("Set thr=3 → variants present in 3–4 samples only (≥5 variants excluded)")
                ], style={'marginLeft': '18px', 'fontSize': '12px', 'color': '#666'}),
                html.Span(
                    "Each group represents a unique subset of variants (no double-counting).",
                    style={'fontSize': '12px', 'color': '#555'}
                )
            ], style={'marginBottom': '8px'}),

            dbc.Button(
                "➕ Add Group",
                id={'type': 'add-group-btn', 'group_type': 'variant'},
                color="secondary",
                size="sm",
                style={
                    'padding': '2px 6px',
                    'fontSize': '12px',
                    'marginTop': '4px'
                }
            ),
        ], style={
            'backgroundColor': '#f8f9fa',
            'padding': '10px',
            'borderRadius': '6px'
        }),

        # --- 구분선 ---
        html.Hr(style={
            'marginTop': '10px',
            'marginBottom': '10px',
            'borderTop': '1px solid #ddd'
        }),

        # --- 그룹 리스트 ---
        html.Div([
            create_group_block(11, visible=True, removable=False, group_type="variant"),
            create_group_block(12, visible=True, removable=False, group_type="variant"),
            create_group_block(13, visible=False, group_type="variant"),
            create_group_block(14, visible=False, group_type="variant"),
            create_group_block(15, visible=False, group_type="variant"),
        ], style={
            'display': 'flex',
            'flexWrap': 'wrap',
            'gap': '10px',
            'justifyContent': 'flex-start',
            'padding': '5px'
        }),
html.Hr(style={'marginTop': '3px', 'marginBottom': '3px'}),
            html.Div(
                id="group-enable-container",
                style={
                    "display": "flex",
                    "flexDirection": "column",
                    "alignItems": "flex-end",
                    "margin": "10px 15px 5px 10px",
                    "gap": "4px"
                },
                children=[
                    dbc.Button(
                        "Apply Grouping to Plots",
                        id="group-enable-switch",
                        n_clicks=0,
                        color="secondary",
                        size="sm",
                        style={"width": "220px"}
                    )
                ]
            )

    ],
    style={
        'border': '1px solid #ccc',
        'borderRadius': '6px',
        'padding': '12px 15px',
        'backgroundColor': '#ffffff',
        'marginBottom': '10px'
    })

@callback(
    Output({'type': 'group-dropdown', 'group_type': 'variant', 'index': ALL}, 'options'),
    Output({'type': 'group-dropdown', 'group_type': 'variant', 'index': ALL}, 'value'),
    Output("presence-sets-store", "data", allow_duplicate=True),
    Input("filter-meta-store", "data"),
    Input({'type': 'group-dropdown', 'group_type': 'variant', 'index': ALL}, 'value'),
    prevent_initial_call=True
)
def update_variant_group_options(filter_meta, group_values):
    if not filter_meta or "presence_sets" not in filter_meta:
        return [[] for _ in group_values], [None for _ in group_values], None

    normalized_sets = normalize_presence_sets(filter_meta["presence_sets"])
    consumed = set([v for vals in group_values if vals for v in vals])

    options_list, values_out = [], []
    for vals in group_values:
        available = [
            {"label": f"Set thr={s['threshold']} (unique={s['count']})", "value": str(s["threshold"])}
            for s in normalized_sets
            if (str(s["threshold"]) not in consumed or (vals and str(s["threshold"]) in vals))
        ]
        options_list.append(available)

        if vals:
            valid_vals = [v for v in vals if v in [opt["value"] for opt in available]]
            values_out.append(valid_vals if valid_vals else None)
        else:
            values_out.append(None)

    return options_list, values_out, normalized_sets

def normalize_presence_sets(presence_sets):
    """
    presence_sets: threshold별 원래 집합들
    return: 차집합 기반 세트 (서로 disjoint)
    """
    # threshold 높은 순으로 정렬
    sorted_sets = sorted(presence_sets, key=lambda x: x["threshold"], reverse=True)
    
    result_sets = []
    higher_variants = set()
    
    for s in sorted_sets:
        current = set(s["variant_ids"])
        unique_variants = current - higher_variants  # 상위 threshold에서 이미 잡힌 건 제외
        if unique_variants:
            result_sets.append({
                "threshold": s["threshold"],
                "variant_ids": list(unique_variants),
                "total_samples": s["total_samples"],
                "count": len(unique_variants)
            })
        # 상위 집합에 추가
        higher_variants |= current

    return result_sets[::-1]  # 낮은 threshold부터 반환


@callback(
    Output({'type': 'group-dropdown', 'group_type': 'sample', 'index': ALL}, 'options'),
    Output({'type': 'group-dropdown', 'group_type': 'sample', 'index': ALL}, 'value'),
    Output("presence-sets-store", "data", allow_duplicate=True),
    Input("gwas-selected-samples-store", "data"),
    Input({'type': 'group-dropdown', 'group_type': 'sample', 'index': ALL}, 'value'),
    prevent_initial_call=True
)
def update_sample_group_options(samples, group_values):
    if not samples or len(samples) < 3:
        return [[] for _ in group_values], [None for _ in group_values], None

    consumed = set([s for vals in group_values if vals for s in vals])
    options_list, values_out = [], []

    for vals in group_values:
        current = set(vals or [])
        available = [s for s in samples if s not in (consumed - current)]
        opts = [{'label': s, 'value': s} for s in available]
        options_list.append(opts)

        if vals:
            valid_vals = [v for v in vals if v in available]
            values_out.append(valid_vals if valid_vals else None)
        else:
            values_out.append(None)

    return options_list, values_out, None

@callback(
    Output("sample-group-container", "children"),
    Input("gwas-selected-samples-store", "data")
)
def render_sample_grouping(samples):
    if not samples or len(samples) < 3:
        return html.Div(
            "⚠️ Please select at least 3 samples to enable sample grouping.",
            style={'color': 'red', 'fontSize': '13px', 'padding': '6px'}
        )

    return html.Div([
        # --- 상단 설명 + Add 버튼 ---
        html.Div([
            html.Div([
                html.I(className="fas fa-info-circle", style={'color': '#007bff', 'marginRight': '6px'}),
                html.Span("Sample grouping", style={
                    'fontWeight': 'bold',
                    'color': '#007bff',
                    'marginRight': '6px'
                }),
                html.Span(
                    "Requires ≥3 samples to form groups like [a,b] vs [c,d].",
                    style={'fontSize': '12px', 'color': '#555'}
                )
            ], style={'marginBottom': '6px', 'display': 'flex', 'alignItems': 'center'}),

            dbc.Button(
                "➕ Add Group",
                id={'type': 'add-group-btn', 'group_type': 'sample'},
                color="secondary",
                size="sm",
                style={
                    'padding': '2px 6px',
                    'fontSize': '12px',
                    'marginTop': '4px'
                }
            ),
        ], style={
            'backgroundColor': '#f8f9fa',
            'padding': '10px',
            'borderRadius': '6px'
        }),

        # --- 구분선 ---
        html.Hr(style={
            'marginTop': '10px',
            'marginBottom': '10px',
            'borderTop': '1px solid #ddd'
        }),

        # --- 그룹 리스트 ---
        html.Div([
            create_group_block(6, visible=True, removable=False, group_type="sample"),
            create_group_block(7, visible=True, removable=False, group_type="sample"),
            create_group_block(8, visible=False, group_type="sample"),
            create_group_block(9, visible=False, group_type="sample"),
            create_group_block(10, visible=False, group_type="sample"),
        ], style={
            'display': 'flex',
            'flexWrap': 'wrap',
            'gap': '10px',
            'justifyContent': 'flex-start',
            'padding': '5px'
        }),
        html.Hr(style={'marginTop': '3px', 'marginBottom': '3px'}),
            html.Div(
                id="group-enable-container",
                style={
                    "display": "flex",
                    "flexDirection": "column",
                    "alignItems": "flex-end",
                    "margin": "10px 15px 5px 10px",
                    "gap": "4px"
                },
                children=[
                    dbc.Button(
                        "Apply Grouping to Plots",
                        id="group-enable-switch",
                        n_clicks=0,
                        color="secondary",
                        size="sm",
                        style={"width": "220px"}
                    )
                ]
            )
    ],
    style={
        'border': '1px solid #ccc',
        'borderRadius': '6px',
        'padding': '12px 15px',
        'backgroundColor': '#ffffff',
        'marginBottom': '10px'
    })

'''
@callback(
    Output("group-enable-container", "style",allow_duplicate=True),
    Output("group-enable-switch", "value",allow_duplicate=True),
    Input({'type': 'group-dropdown', 'group_type': ALL, 'index': ALL}, "value")
    ,prevent_initial_call=True
)
def toggle_group_enable(all_group_values):
    has_any_value = any(
        v for v in all_group_values if v  # 리스트나 문자열 값이 있을 때 True
    )
    if has_any_value:
        return {"display": "block", "marginBottom": "10px"}, False
    else:
        return {"display": "none"}, False
'''


def create_sample_summary_data_from_combo_bygroup(
    combo_json, vcf_samples,
    maf_cutoff=None, maf_enabled=False,
    presence_threshold=None, presence_enabled=False,
    groupset=None,
    presence_sets=None
):
    import pandas as pd

    if not combo_json or combo_json.get("status") != "ready":
        return None, None, {}

    df_vo = _records_to_df(combo_json.get("variant_only"))
    df_vp = _records_to_df(combo_json.get("variant_pvalue"))
    if df_vo.empty:
        return None, df_vp, {}

    if not df_vp.empty and "Subtrait" in df_vp.columns:
        df_vp["Clean_Subtrait"] = df_vp["Subtrait"].apply(remove_parentheses_numbers)

    # --- MAF filter ---
    if maf_enabled and maf_cutoff is not None:
        df_vp = df_vp[df_vp["MAF"].astype(float) >= maf_cutoff]
        valid_variants = set(df_vp["Variant ID"].unique())
        df_vo = df_vo[df_vo["Variant ID"].isin(valid_variants)]

    rows = []
    df_vp_grouped = []
    mode = groupset.get("mode") if groupset else None
    gvals_all = {str(k): v for k, v in (groupset.get("values") or {}).items()}

    # ✅ 그룹 인덱스 구간 정의
    mode_ranges = {
        "trait": range(1, 6),
        "sample": range(6, 11),
        "variant": range(11, 16)
    }

    # ✅ 현재 모드의 그룹만 필터링
    target_range = mode_ranges.get(mode, range(1, 6))
    gvals = {k: v for k, v in gvals_all.items()
             if int(k.replace("group", "")) in target_range}

    # ✅ normalize_group_label (1~5 라벨로 환원)
    def normalize_group_label(gname):
        try:
            gnum = int(gname.replace("group", ""))
        except:
            return gname
        if mode == "trait":
            newnum = gnum
        elif mode == "sample":
            newnum = gnum - 5
        elif mode == "variant":
            newnum = gnum - 10
        else:
            newnum = gnum
        return f"group{newnum}"

    # ---------------------------
    # (이 아래부터 기존 로직 동일)
    # ---------------------------

    # ----------------- Variant mode -----------------
    if mode == "variant" and presence_sets:
        set_map = {str(s["threshold"]): set(s["variant_ids"]) for s in presence_sets}
        for gname, thresholds in gvals.items():
            if not thresholds:
                continue
            variant_ids = set().union(*(set_map[str(t)] for t in thresholds if str(t) in set_map))
            group_label = normalize_group_label(gname)

            sdf = df_vo[df_vo["Variant ID"].isin(variant_ids)]
            for trait, tdf in sdf.groupby("Trait"):
                rows.append({
                    "Trait": trait,
                    "Group": group_label,
                    "Variant Count": int(tdf["Variant ID"].nunique())
                })

            if not df_vp.empty:
                df_vp_masked = df_vp[df_vp["Variant ID"].isin(variant_ids)].copy()
                df_vp_masked["Group"] = group_label
                df_vp_grouped.append(df_vp_masked)
        df_vp_final = pd.concat(df_vp_grouped, ignore_index=True) if df_vp_grouped else pd.DataFrame()

    # ----------------- Sample mode -----------------
    elif mode == "sample":
        norm_gvals = {normalize_group_label(g): (vals or []) for g, vals in gvals.items()}
        participating_samples = sorted({s for lst in norm_gvals.values() for s in lst})

        col_of = lambda s: f"{s}_GT"
        all_part_gt_cols = [col_of(s) for s in participating_samples if col_of(s) in df_vo.columns]

        group_sets = {}

        for gname, sample_list in norm_gvals.items():
            if not sample_list:
                continue
            group_cols = [col_of(s) for s in sample_list if col_of(s) in all_part_gt_cols]
            if not group_cols:
                continue
            outside_samples = sorted({s for og, lst in norm_gvals.items() if og != gname for s in lst})
            outside_cols = [col_of(s) for s in outside_samples if col_of(s) in all_part_gt_cols]

            mask_all_in_group = df_vo[group_cols].notna().all(axis=1)
            mask_none_outside = ~df_vo[outside_cols].notna().any(axis=1) if outside_cols else True

            vids = set(df_vo.loc[mask_all_in_group & mask_none_outside, "Variant ID"].unique())
            group_sets[gname] = vids

        df_vp_grouped = []
        for gname, vids in group_sets.items():
            if not vids:
                continue
            group_label = f"{gname}_only"

            sdf = df_vo[df_vo["Variant ID"].isin(vids)]
            for trait, tdf in sdf.groupby("Trait"):
                rows.append({
                    "Trait": trait,
                    "Group": group_label,
                    "Variant Count": int(tdf["Variant ID"].nunique())
                })

            if not df_vp.empty:
                df_vp_masked = df_vp[df_vp["Variant ID"].isin(vids)].copy()
                df_vp_masked["Group"] = group_label
                df_vp_grouped.append(df_vp_masked)

        df_vp_final = pd.concat(df_vp_grouped, ignore_index=True) if df_vp_grouped else pd.DataFrame()

    # ----------------- Trait mode -----------------
    elif mode == "trait":
        for gname, traits in gvals.items():
            if not traits:
                continue
            group_label = normalize_group_label(gname)
            tdf = df_vo[df_vo["Trait"].isin(traits)]
            for trait, sdf in tdf.groupby("Trait"):
                rows.append({
                    "Trait": trait,
                    "Group": group_label,
                    "Variant Count": int(sdf["Variant ID"].nunique())
                })
            if not df_vp.empty:
                df_vp_masked = df_vp[df_vp["Trait"].isin(traits)].copy()
                df_vp_masked["Group"] = group_label
                df_vp_grouped.append(df_vp_masked)
        df_vp_final = pd.concat(df_vp_grouped, ignore_index=True) if df_vp_grouped else pd.DataFrame()
    else:
        df_vp_final = pd.DataFrame()

    df_summary = pd.DataFrame(rows)
    return df_summary, df_vp_final, {}

def extract_group_index(group_name):
    import re
    match = re.search(r"group[_\-]?(\d+)", group_name)
    return int(match.group(1)) if match else float('inf')


def create_enhanced_clickable_barplot_bygroup(summary_df):
    """
    그룹 모드용 barplot
    - trait 모드: 그룹별 색만 다르게 표시 (grouped bar)
    - variant / sample 모드: stacked bar 형식 + group 색
    """
    from dash import html, dcc
    import plotly.express as px

    if summary_df is None or summary_df.empty:
        return html.Div("No group summary data available", className="text-muted text-center p=3")

    # mode 추론 (Trait 모드인지, 나머지인지)
    mode = "trait"
    if "Group" in summary_df.columns:
        # heuristic: 그룹 이름에 "_only" 있거나 threshold 기반이면 variant/sample 모드
        unique_groups = summary_df["Group"].unique()
        if any("_only" in g for g in unique_groups) or any("Set" in str(g) or "thr" in str(g) for g in unique_groups):
            mode = "stacked"
        else:
            mode = "trait"

    # 색상 팔레트
    color_discrete_map = {}
    palette = px.colors.qualitative.Set2
    sorted_groups = sorted(summary_df["Group"].unique(), key=extract_group_index)

    for i, g in enumerate(sorted_groups):
        color_discrete_map[g] = palette[i % len(palette)]

    # px.bar 생성
    fig = px.bar(
        summary_df,
        x="Trait",
        y="Variant Count",
        color="Group",
        text="Variant Count",
        color_discrete_map=color_discrete_map
    )

    if mode == "trait":
        fig.update_layout(barmode="group")
    else:
        fig.update_layout(barmode="stack")

    # 공통 레이아웃
    fig.update_traces(textposition='auto')
    fig.update_layout(
        title="Trait Barplot (Group Mode)",
        xaxis_title="Traits",
        yaxis_title="Variant Count",
        xaxis_tickangle=45,
        height=500,
        margin=dict(l=50, r=30, t=60, b=100),
        plot_bgcolor='white',
        paper_bgcolor='white',
        clickmode='event',
        hovermode='closest',
    )
    #fig.update_xaxes(showgrid=True, gridcolor='#f0f0f0')
    fig.update_xaxes(visible=False, showticklabels=False, showgrid=False)
    fig.update_yaxes(showgrid=True, gridcolor='#f0f0f0')

    return fig




@callback(
    [
        Output("integrated-trait-barplot-group", "figure"),
        Output("filter-meta-store2", "data"),
        Output("df-vp-store2", "data"),
        Output("group-info-display", "children"),
        Output("group-info-display", "style"),
    ],
    [
        Input("gwas-selected-samples-store", "data"),
        Input("integrated-subtrait-dropdown", "value"),
        Input("gwas-combo-store", "data"),
        # 🔽 maf-filter-store → State로 이동
        Input("snp-occurrence-store", "data"),
        Input("group-enable-switch", "n_clicks"),
        Input("group-type-radio", "value"),
        Input({'type': 'group-dropdown', 'group_type': ALL, 'index': ALL}, 'value'),
        # 🔽 새로 maf-enabled만 Input으로 추가 (on/off 반응)
        Input("maf-enabled", "value"),
    ],
    [
        State("maf-filter-store", "data"),
        State("filter-meta-store2", "data"),
        State("presence-sets-store", "data"),
        State({'type': 'group-dropdown', 'group_type': ALL, 'index': ALL}, 'id')
    ],
    prevent_initial_call=True
)
def update_trait_bar_container_group(
    selected_samples, selected_subtraits, combo_json,
    presence_state, group_enabled, group_mode, group_values,
    maf_enabled,    # 새 input
    maf_filter_state, filter_meta_state, presence_sets_data, group_ids
):
    print("🔄 Updating group trait bar container...")

    if not group_enabled:
        return dash.no_update, dash.no_update, dash.no_update,dash.no_update,dash.no_update

    selected_samples = selected_samples or []
    selected_subtraits = selected_subtraits or []
    filter_meta_state = filter_meta_state or {}

    # --- 기본 확인 ---
    if not selected_samples:
        return html.Div("Please select samples", style={'textAlign': 'center', 'padding': '20px'}), filter_meta_state, None, dash.no_update,dash.no_update
    if not GWAS_AVAILABLE:
        return html.Div("GWAS module not available", style={'textAlign': 'center', 'padding': '20px'}), filter_meta_state, None, dash.no_update,dash.no_update
    if not combo_json or combo_json.get("status") != "ready":
        return html.Div("GWAS data not ready", style={'textAlign': 'center', 'padding': '20px'}), filter_meta_state, None,dash.no_update,dash.no_update

    # --- 필터 상태 ---
    maf_cutoff_default = 0.05
    maf_cutoff = maf_cutoff_default
    maf_active = False

    if maf_enabled:  # ✅ 스위치 켜진 경우만 활성화
        if maf_filter_state:
            maf_cutoff = maf_filter_state.get("cutoff", maf_cutoff_default)
            maf_active = bool(maf_filter_state.get("enabled", True))
        else:
            maf_cutoff = maf_cutoff_default
            maf_active = True
    else:
        maf_active = False
        maf_cutoff = maf_cutoff_default 



    presence_enabled = presence_state.get("enabled", False) if presence_state else False
    presence_threshold = presence_state.get("threshold", None) if presence_state else None

    try:
        # --- 그룹셋 구성 ---
        groupset_values = {}
        for id_obj, vals in zip(group_ids, group_values):
            # id_obj 예: {'type': 'group-dropdown', 'group_type': 'sample', 'index': 6}
            gtype = id_obj.get('group_type')
            gidx = id_obj.get('index')
            groupset_values[f"group{gidx}"] = vals

        groupset = {
            "mode": group_mode,
            "values": groupset_values
        }
        # --- 그룹 기반 데이터 생성 ---
        summary_df, df_vp2, filter_meta_new = create_sample_summary_data_from_combo_bygroup(
            combo_json, selected_samples,
            maf_cutoff=maf_cutoff, maf_enabled=maf_enabled,
            presence_threshold=presence_threshold, presence_enabled=presence_enabled,
            groupset=groupset,
            presence_sets=presence_sets_data or []
        )

        if summary_df is None or summary_df.empty:
            return html.Div("No group summary data", style={'textAlign': 'center', 'padding': '20px'}), filter_meta_state, None, dash.no_update,dash.no_update

        # --- Subtrait 필터 ---
        df_plot = summary_df
        df_vp = df_vp2
        if selected_subtraits and "Subtrait" in summary_df.columns:
            df_plot = summary_df[summary_df["Subtrait"].isin(selected_subtraits)].copy()
            df_vp = df_vp2[df_vp2["Subtrait"].isin(selected_subtraits)].copy()
            if df_plot.empty:
                return html.Div(
                    f"No traits for subtraits: {', '.join(selected_subtraits)}",
                    style={'textAlign': 'center', 'padding': '20px'}
                ), filter_meta_state, None, dash.no_update,dash.no_update

        # --- Presence set 누적 저장 ---
        MODE_GROUP_RANGES = {
            "trait": range(1, 6),      # group1 ~ group5
            "sample": range(6, 11),    # group6 ~ group10
            "variant": range(11, 16)   # group11 ~ group15
        }
        display_children = []

        for idx, (gname, vals) in enumerate(groupset["values"].items(), start=1):
            if not vals:
                continue

            # 현재 모드에 맞는 그룹 범위가 아니면 skip
            if idx not in MODE_GROUP_RANGES[group_mode]:
                continue

            # ─── Trait 모드 ───
            if group_mode == "trait":
                label = format_group_label("trait", idx)
                display_children.append(
                    html.Div(f"{label}: {', '.join(vals)}",
                            style={"fontSize": "12px", "marginBottom": "4px"})
                )

            # ─── Sample 모드 ───
            elif group_mode == "sample":
                sample_label = format_group_label("sample", idx - 5)  # 6→1, 7→2 ...
                display_children.append(
                    html.Div(f"{sample_label}: {', '.join(vals)}",
                            style={"fontSize": "12px", "marginBottom": "4px"})
                )

            # ─── Variant 모드 ───
            elif group_mode == "variant":
                labels = []
                for v in vals:
                    label = next(
                        (f"Set thr={s['threshold']} (unique={s['count']})"
                        for s in presence_sets_data
                        if str(s["threshold"]) == v),
                        f"thr={v}"
                    )
                    labels.append(label)

                group_label = format_group_label("variant", idx - 10)  # 11→1, 12→2 ...
                display_children.append(
                    html.Div(f"{group_label}: {', '.join(labels)}",
                            style={"fontSize": "12px", "marginBottom": "4px"})
                )

        display_style = {
                        "fontSize": "12px",
                        "color": "#444",
                        "margin": "10px auto",
                        "padding": "6px",
                        "border": "1px solid #ddd",
                        "borderRadius": "4px",
                        "width": "90%",
                        "display": "block"} if display_children else {"display": "none"}

        if presence_enabled and filter_meta_new and "presence" in filter_meta_new:
            presence_info = filter_meta_new["presence"]
            presence_sets = filter_meta_state.get("presence_sets", [])
            if not any(p["threshold"] == presence_info["threshold"] for p in presence_sets):
                presence_sets.append(presence_info)
            filter_meta_state["presence_sets"] = presence_sets

        # --- Barplot 생성 ---
        df_plot = df_plot.drop(columns=["Subtrait"], errors="ignore")
        df_plot = df_plot.rename(columns={"Clean_Subtrait": "Subtrait"})
        graph = create_enhanced_clickable_barplot_bygroup(df_plot)

        return graph, filter_meta_state, df_vp.to_dict("records"), display_children,display_style

    except Exception as e:
        print(f"❌ Error in group bar container callback: {e}")
        return html.Div(f"Error: {e}", style={'color': 'red', 'textAlign': 'center'}), filter_meta_state, None, dash.no_update,dash.no_update

@callback(
    #Output('group-enable-switch', 'value', allow_duplicate=True),
    Output('subtrait-container-wrapper-filter', 'style'),
    Output('subtrait-container-wrapper-group', 'style'),
    Output('close-group-barplot-btn', 'style'),   # ✅ 버튼 표시 제어 추가
    Input('integrated-trait-barplot-group', 'figure'),
    prevent_initial_call=True
)
def toggle_barplot_display_by_group(fig):
    """
    그룹 figure 존재 여부에 따라 Filter / Group plot 표시 전환
    + 그룹 plot 시 Close 버튼 노출
    """
    base_btn_style = {
        #"position": "absolute",
        "top": "8px",
        "right": "12px",
        "border": "none",
        "background": "transparent",
        "color": "#666",
        "fontSize": "16px",
        "cursor": "pointer",
        #"zIndex": 5,
    }

    # 🔹 그룹 플롯이 없는 경우 → 기본 barplot 모드
    if not fig or not fig.get("data"):
        return (
            #False,  # group-enable-switch OFF
            {'display': 'block'},  # filter 영역 보임
            {'display': 'none'},   # group 영역 숨김
            {**base_btn_style, "display": "none"}  # ❌ 버튼 숨김
        )

    # 🔹 그룹 플롯 존재 → group plot 활성
    return (
        #False,
        {'display': 'none'},       # filter 숨김
        {'display': 'block'},      # group 영역 표시
        {**base_btn_style, "display": "block"}  # ✅ 버튼 표시
    )

@callback(
    Output('subtrait-container-wrapper-filter', 'style', allow_duplicate=True),
    Output('subtrait-container-wrapper-group', 'style', allow_duplicate=True),
    Output('gwas-sample-scatter-container', 'style', allow_duplicate=True),
    Output('gwas-sample-scatter-group-container', 'style', allow_duplicate=True),
    Output('close-group-barplot-btn', 'style', allow_duplicate=True),
    Input('close-group-barplot-btn', 'n_clicks'),
    prevent_initial_call=True
)
def close_group_barplot_all(n_clicks):
    """
    ❌ 버튼 클릭 시:
    - Group plot 종료 → Filter plot 복귀
    - subtrait 및 gwas 영역 모두 filter 모드로 되돌림
    """
    from dash import exceptions

    if not n_clicks:
        raise exceptions.PreventUpdate

    # 공통 버튼 스타일
    base_btn_style = {
        "border": "none",
        "background": "transparent",
        "color": "#666",
        "fontSize": "14px",
        "cursor": "pointer",
        "marginTop": "4px",
    }

    return (
        {'display': 'block'},  # subtrait: filter 보이기
        {'display': 'none'},   # subtrait: group 숨기기
        {'display': 'block'},  # gwas: filter 탭 보이기
        {'display': 'none'},   # gwas: group 탭 숨기기
        {**base_btn_style, "display": "none"}  # 닫기 버튼 숨김
    )


@dash.callback(
    Output('gwas-sample-scatter-container', 'style'),
    Output('gwas-sample-scatter-group-container', 'style'),
    Input('integrated-trait-barplot-group', 'figure'),
    prevent_initial_call=True
)
def toggle_gwas_tab_by_group(fig):
    """
    그룹 Figure 유무에 따라
    - 일반탭(gwas-sample-scatter-container)
    - 그룹탭(gwas-sample-scatter-group-container)
      표시를 전환한다.
    """

    # 기본 숨김값
    filter_style = {'display': 'none'}
    group_style = {'display': 'none'}

    # 🔹 그룹 플롯이 존재하지 않으면 → 기본(Filter-based) 탭 보이기
    if not fig or not fig.get("data"):
        filter_style = {'display': 'block'}
        group_style = {'display': 'none'}
    else:
        # 🔹 그룹 플롯 존재 → Group 탭 보이기
        filter_style = {'display': 'none'}
        group_style = {'display': 'block'}

    return filter_style, group_style



def convert_unique_to_wide_format(unique_dict):
    import pandas as pd

    records = []
    for sample, rows in unique_dict.items():
        for row in rows:
            row = row.copy()
            row[f"{sample}_GT"] = row.get(f"{sample}_GT", row.get("GT", ""))
            row["Sample"] = sample
            records.append(row)

    df_long = pd.DataFrame(records)

    # 공통 key
    base_cols = ['Variant ID','Chromosome','Position','Minor Allele','MAF',
                 'P-value','PMID','Trait',"Subtrait",'Clean_Subtrait']
    
    gt_cols = [col for col in df_long.columns if col.endswith("_GT")]

    df_pivot = df_long.pivot_table(
        index=base_cols,
        values=gt_cols,
        aggfunc="first"
    ).reset_index()

    return df_pivot

@callback(
    Output("gwas-sample-table", "data"),
    Output("gwas-sample-table", "columns"),
    Output("gwas-table-status", "children"),  # ✅ 상태 표시용
    Input("final-df-store", "data"),
    State("sample-mode-store", "data"),
    prevent_initial_call=True
)
def update_table(final_df, mode):
    
    import pandas as pd

    if not final_df:
        return [], [], "⚠ No data available."

    # ✅ UNIQUE 모드 제한
    if mode == "unique":
        n_samples = len(final_df.keys()) if isinstance(final_df, dict) else 0
        total_records = sum(len(rows) for rows in final_df.values())
        if n_samples > 4:
            return [], [], "❌ UNIQUE mode supports ≤4 samples only. Table is not available."
        if total_records == 0:
                return [], [], "⚠ No unique variants available."
                
    if mode == "unique":
        df = convert_unique_to_wide_format(final_df)
    else:
        df = pd.DataFrame(final_df)

    if df.empty:
        return [], [], "⚠ Empty dataset."

    # ✅ 컬럼 이름 처리
    if "Clean_Subtrait" in df.columns:
        df = df.drop(columns=["Subtrait"], errors="ignore")
        df = df.rename(columns={"Clean_Subtrait": "Subtrait"})

    # ✅ 정렬
    df["Chromosome"] = pd.to_numeric(df["Chromosome"], errors="coerce")
    df["Position"] = pd.to_numeric(df["Position"], errors="coerce")
    df = df.sort_values(by=["Chromosome", "Position"])

    columns = [{"name": col, "id": col} for col in df.columns]
    return df.to_dict("records"), columns, ""


@callback(
    [Output("gwas-table-download", "data"),
     Output("gwas-table-download-status", "children")],   # ✅ 상태 메세지 출력용
    Input("gwas-table-download-btn", "n_clicks"),
    State("final-df-store", "data"),
    State("sample-mode-store", "data"),
    prevent_initial_call=True
)
def download_table(n_clicks, final_df, mode):
    import pandas as pd

    if not final_df:
        return dash.no_update, "⚠ No data available."

    # ✅ UNIQUE 모드 제한
    if mode == "unique":
        n_samples = len(final_df.keys()) if isinstance(final_df, dict) else 0
        total_records = sum(len(rows) for rows in final_df.values())
        if n_samples > 4:
            return dash.no_update, "❌ UNIQUE mode supports ≤4 samples only. Table download is not available."
        if total_records == 0:
            return dash.no_update, "⚠ No unique variants available."
        
    # ✅ 데이터 변환
    if mode == "unique":
        df = convert_unique_to_wide_format(final_df)
    else:
        df = pd.DataFrame(final_df)

    # ✅ 컬럼 정리
    if "Clean_Subtrait" in df.columns:
        df = df.drop(columns=["Subtrait"], errors="ignore")
        df = df.rename(columns={"Clean_Subtrait": "Subtrait"})

    # ✅ 정렬
    df["Chromosome"] = pd.to_numeric(df["Chromosome"], errors="coerce")
    df["Position"] = pd.to_numeric(df["Position"], errors="coerce")
    df = df.sort_values(by=["Chromosome", "Position"])

    return dcc.send_data_frame(df.to_csv, "gwas_table.csv", index=False), ""


@callback(
    Output("integrated-filter-section", "style"),
    Output("filter-collapse-btn", "children"),
    Input("filter-collapse-btn", "n_clicks"),
    State("integrated-filter-section", "style"),
    prevent_initial_call=True
)
def toggle_integrated_filter(n, current_style):
    if not n:
        raise PreventUpdate

    is_open = current_style.get("display", "block") == "block"

    if is_open:
        # 닫기
        return (
            {"display": "none", "margin-bottom": "15px", "max-width": "100%"},
            [html.I(className="fas fa-chevron-down", style={"marginRight": "6px"}), "Show Filters"]
        )
    else:
        # 다시 펼치기
        return (
            {"display": "block", "margin-bottom": "15px", "max-width": "100%"},
            [html.I(className="fas fa-chevron-up", style={"marginRight": "6px"}), "Hide Filters"]
        )

@callback(
    Output('gwas-sample-scatter-group', 'figure'),
    Output('gwas-sample-scatter-group-data', 'data'),  # 🆕 used_samples 저장용
    Output('gwas-sample-scatter-group-container', 'style',allow_duplicate=True),
    
    Input('df-vp-store2', 'data'),
    Input('gwas-click-store', 'data'),
    #Input('scatter-refresh-trigger-group', 'data')]
    
    prevent_initial_call=True
)
def update_group_scatter(df_vp_data, click_store):
    import pandas as pd, numpy as np
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    from dash import ctx
    trig = ctx.triggered_id



    if not df_vp_data:
        return go.Figure().add_annotation(
            text="No group data available",
            x=0.5, y=0.5, xref="paper", yref="paper",
            showarrow=False, font=dict(size=16, color="gray")
        ), no_update, {'display': 'none'}

    df = pd.DataFrame(df_vp_data)
    if df.empty or "Group" not in df.columns:
        return go.Figure().add_annotation(
            text="No valid group data",
            x=0.5, y=0.5, xref="paper", yref="paper",
            showarrow=False, font=dict(size=16, color="gray")
        ), no_update, {'display': 'none'}
    
    if trig == 'gwas-click-store':
    # None -> 패스, []만 재생성, 길이>0 -> 패스
        if click_store is None or len(click_store) != 0:
            raise PreventUpdate

    def format_gt_row(row):
        if not gt_cols:
            return None
        if len(gt_cols) == 1:
            return row[gt_cols[0]]
        # 여러개인 경우 "Sample1:GT1, Sample2:GT2" 문자열로 변환
        parts = [f"{c.replace('_GT','')}: {row[c]}" for c in gt_cols if pd.notna(row[c])]
        return ", ".join(parts)
    
    gt_cols = [c for c in df.columns if c.endswith("_GT")]
    df["GT_Info"] = df.apply(format_gt_row, axis=1)
    gt_label = "GT" if len(gt_cols) == 1 else "GTs"
    df["Group_Info"] = df["Group"] if "Group" in df.columns else ""






    df["P-value"] = pd.to_numeric(df["P-value"], errors="coerce")
    df["neg_log10_pvalue"] = df["P-value"].apply(lambda x: -np.log10(x) if pd.notna(x) and x > 0 else None)

    chromosomes = sorted(df["Chromosome"].dropna().unique(), key=lambda x: int(x))
    cols = min(4, len(chromosomes))
    rows = (len(chromosomes) + cols - 1) // cols
    subplot_titles = [f"CHR{c}" for c in chromosomes]

    fig = make_subplots(
        rows=rows, cols=cols,
        subplot_titles=subplot_titles,
        vertical_spacing=0.12,
        horizontal_spacing=0.08
    )

    color_palette = px.colors.qualitative.Set2
    sorted_groups = sorted(df["Group"].unique())
    color_map = {group: color_palette[i % len(color_palette)] for i, group in enumerate(sorted_groups)}

    seen_legends = set()
    for i, chr_id in enumerate(chromosomes):
        row, col = i // cols + 1, i % cols + 1
        subdf = df[df["Chromosome"] == chr_id]
        if subdf.empty:
            continue

        for group, grp in subdf.groupby("Group"):
            showlegend = group not in seen_legends
            fig.add_trace(go.Scattergl(
                x=grp["Position"], y=grp["neg_log10_pvalue"],
                mode="markers",
                marker=dict(color=color_map[group], size=7, line=dict(width=1, color="black")),
                name=group,
                legendgroup=group,
                showlegend=showlegend,
                customdata=np.stack([grp["Chromosome"], 
                grp["P-value"], 
                grp["Subtrait"],
                grp["Minor Allele"],
                grp["Variant ID"],
                grp["GT_Info"],
                grp["Group_Info"]
                ],  axis=-1),
                hovertemplate=(
                    "Chr: %{customdata[0]}<br>"
                    "Pos: %{x}<br>"
                    "P: %{customdata[1]:.2e}<br>"
                    "-log₁₀(P): %{y:.2f}<br>"
                    "Subtrait: %{customdata[2]}<br>"
                    "Minor Allele: %{customdata[3]}<br>"
                    "Variant ID: %{customdata[4]}<br>"
                    + gt_label + ": %{customdata[5]}<br>"
                    "Group: %{customdata[6]}<extra></extra>"
                )
            ), row=row, col=col)
            seen_legends.add(group)

    fig.update_layout(
        height=800,
        showlegend=True,
        margin=dict(t=50, b=40, l=50, r=20),
        xaxis_title="Position",
        yaxis_title="-log₁₀(P)"
    )

    used_samples = [c.replace("_GT", "") for c in gt_cols]
    return fig, {"used_samples": used_samples},  {'display': 'block'}

@callback(
    Output("gwas-sample-table-group", "data"),
    Output("gwas-sample-table-group", "columns"),
    Input("df-vp-store2", "data"),
    prevent_initial_call=True
)
def update_group_table(data):
    import pandas as pd
    if not data:
        return [], []

    df = pd.DataFrame(data)
    if df.empty:
        return [], []
    
    if "Clean_Subtrait" in df.columns:
        df = df.drop(columns=["Subtrait"], errors="ignore")
        df = df.rename(columns={"Clean_Subtrait": "Subtrait"})
    
    df = df.drop(columns=["Group_Info", "GT_Info"], errors="ignore")

    df["Chromosome"] = pd.to_numeric(df["Chromosome"], errors="coerce")
    df["Position"] = pd.to_numeric(df["Position"], errors="coerce")
    df = df.sort_values(by=["Chromosome", "Position"])

    
    columns = [{"name": col, "id": col} for col in df.columns]
    return df.to_dict("records"), columns

@callback(
    Output("gwas-table-download-group", "data"),
    Input("gwas-table-download-btn-group", "n_clicks"),
    State("df-vp-store2", "data"),
    prevent_initial_call=True
)
def download_group_table(n_clicks, data):
    
    print("download_group_table",n_clicks)
    
    if not data:
        return dash.no_update
    
    df = pd.DataFrame(data)
    if df.empty:
        return [], []

    if "Clean_Subtrait" in df.columns:
        df = df.drop(columns=["Subtrait"], errors="ignore")
        df = df.rename(columns={"Clean_Subtrait": "Subtrait"})
    
    df = df.drop(columns=["Group_Info", "GT_Info"], errors="ignore")
    df["Chromosome"] = pd.to_numeric(df["Chromosome"], errors="coerce")
    df["Position"] = pd.to_numeric(df["Position"], errors="coerce")
    df = df.sort_values(by=["Chromosome", "Position"])
    return dcc.send_data_frame(df.to_csv, "group_gwas_table.csv", index=False)

'''
@callback(
    Output('gwas-sample-scatter-wrapper', 'style',allow_duplicate=True),
    Output('gwas-sample-scatter-group-wrapper', 'style',allow_duplicate=True),
    Output('gwas-sample-table-wrapper', 'style',allow_duplicate=True),
    Output('gwas-sample-table-group-wrapper', 'style'),
    Input('gwas-plot-mode-radio', 'value'),
    prevent_initial_call=True
)
def toggle_gwas_group_filter_views(mode):
    """
    라디오 선택에 따라 적절한 scatter/table container 보여줌.
    - scatter: 기본 filter 기반 scatter
    - table: filter 기반 table
    - Gscatter: group 기반 scatter
    - Gtable: group 기반 table
    """
    # 기본 전부 숨김
    f_scatter = {'display': 'none'}
    f_table = {'display': 'none'}
    g_scatter = {'display': 'none'}
    g_table = {'display': 'none'}

    if mode == "scatter":
        f_scatter = {'display': 'block'}
    elif mode == "table":
        f_table = {'display': 'block'}
    elif mode == "Gscatter":
        g_scatter = {'display': 'block'}
    elif mode == "Gtable":
        g_table = {'display': 'block'}

    return f_scatter, g_scatter, f_table, g_table
'''


#==============================================================

_ID_COL = "id"

def get_category_color_theme(category_name):
    """카테고리별 일관된 색상 테마 반환 (원형 유지)"""
    category_themes = {
        'yield_quality': {
            'primary': '#e74c3c',
            'secondary': '#c0392b',
            'accent': '#ec7063',
            'light': '#fadbd8',
            'gradient_start': '#e74c3c',
            'gradient_end': '#f1948a'
        },
        'cultivation': {
            'primary': '#f39c12',
            'secondary': '#d68910',
            'accent': '#f8c471',
            'light': '#fdebd0',
            'gradient_start': '#f39c12',
            'gradient_end': '#f7dc6f'
        },
        'disaster_resistance': {
            'primary': '#27ae60',
            'secondary': '#229954',
            'accent': '#58d68d',
            'light': '#d5f4e6',
            'gradient_start': '#27ae60',
            'gradient_end': '#82e5aa'
        },
        'morphological': {
            'primary': '#3498db',
            'secondary': '#2980b9',
            'accent': '#5dade2',
            'light': '#d6eaf8',
            'gradient_start': '#3498db',
            'gradient_end': '#85c1e9'
        }
    }
    return category_themes.get(category_name, {
        'primary': '#6c757d',
        'secondary': '#495057',
        'accent': '#adb5bd',
        'light': '#f8f9fa',
        'gradient_start': '#6c757d',
        'gradient_end': '#adb5bd'
    })

def get_enhanced_seaborn_palette(category_name, num_colors=15):
    """카테고리별 최적화된 색상 팔레트 (원형 유지)"""
    if category_name == 'yield_quality':
        return [
            '#e74c3c', '#f39c12', '#f1c40f', '#e67e22', '#d35400',
            '#c0392b', '#a93226', '#922b21', '#7b241c', '#641e16',
            '#ff7f50', '#ff6347', '#ff4500', '#ffa500', '#ffd700'
        ][:num_colors]
    elif category_name == 'cultivation':
        return [
            '#f39c12', '#27ae60', '#16a085', '#d4ac0d', '#b7950b',
            '#a04000', '#935116', '#7e5109', '#6e2c00', '#5d4e75',
            '#daa520', '#cd853f', '#d2691e', '#bc8f8f', '#f4a460'
        ][:num_colors]
    elif category_name == 'disaster_resistance':
        return [
            '#27ae60', '#2980b9', '#8e44ad', '#16a085', '#2e8b57',
            '#229954', '#1f4e79', '#76448a', '#138d75', '#0e6b0e',
            '#006400', '#483d8b', '#4682b4', '#008b8b', '#556b2f'
        ][:num_colors]
    elif category_name == 'morphological':
        return [
            '#3498db', '#9b59b6', '#1abc9c', '#34495e', '#5d6d7e',
            '#2980b9', '#8e44ad', '#16a085', '#2c3e50', '#566573',
            '#4169e1', '#6a5acd', '#20b2aa', '#708090', '#778899'
        ][:num_colors]
    else:
        return [
            '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
            '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
            '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5'
        ][:num_colors]

def detect_trait_category(trait_name):
    """트레이트 이름으로부터 카테고리 자동 감지 (원형 유지)"""
    trait_lower = (trait_name or "").lower()
    yield_quality_keywords = [
        'endosperm', 'hull', 'grain', 'weight', 'amylose', 'protein',
        'alkali', 'length', 'width', 'shape', 'core', 'belly', 'panicle',
        'phenol', 'anthocyanin', 'scavenging', 'yield', 'quality'
    ]
    cultivation_keywords = [
        'date', 'planting', 'seeding', 'heading', 'variety_type',
        'season', 'time', 'panicles_per_plant'
    ]
    resistance_keywords = [
        'resistance', 'disease', 'blight', 'blast', 'stress', 'cold',
        'bacterial', 'stripe', 'brown', 'planthopper', 'pit', 'temp',
        'germination', 'viviparous', 'bakanae', 'shattering', 'pi'
    ]
    morphological_keywords = [
        'color', 'tillering', 'habit', 'internode', 'flag', 'angle',
        'exertion', 'culm', 'length', 'awn', 'sheath', 'erectness',
        'apiculus', 'leaf', 'morphology'
    ]
    scores = {
        'yield_quality': sum(kw in trait_lower for kw in yield_quality_keywords),
        'cultivation': sum(kw in trait_lower for kw in cultivation_keywords),
        'disaster_resistance': sum(kw in trait_lower for kw in resistance_keywords),
        'morphological': sum(kw in trait_lower for kw in morphological_keywords)
    }
    detected = max(scores.items(), key=lambda x: x[1])[0]
    return detected if scores[detected] > 0 else 'morphological'

def _coerce_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.astype(str).str.strip(), errors='coerce')

def get_trait_type_from_db(trait_name):
    """
    ✅ DB 대신 PHENO_DF 기반으로 타입 추론
    - 80% 이상 숫자 파싱 성공 → numeric
    - 컬럼명이 date/mmdd 포함 → date
    - 그 외 → categorical
    """
    if PHENO_DF is None or PHENO_DF.empty or trait_name not in PHENO_DF.columns:
        return "categorical"
    s = PHENO_DF[trait_name].dropna().astype(str)
    numeric_ratio = _coerce_numeric(s).notna().mean()
    if ("date" in trait_name.lower()) or ("mmdd" in trait_name.lower()):
        return "date"
    if numeric_ratio >= 0.8:
        return "numeric"
    return "categorical"

def get_trait_full_range_data(trait_name):
    """
    ✅ DB 대신 PHENO_DF로 전체 범위/통계 계산 (Global Scale 보장)
    반환: {min,max,count,avg,stddev} 또는 None
    """
    if PHENO_DF is None or PHENO_DF.empty or trait_name not in PHENO_DF.columns:
        return None
    s = _coerce_numeric(PHENO_DF[trait_name]).dropna()
    if s.empty:
        return None
    return {
        "min": float(s.min()),
        "max": float(s.max()),
        "count": int(s.size),
        "avg": float(s.mean()),
        "stddev": float(s.std())
    }

def generate_enhanced_gradient_color(norm_value, theme, category_name):
    """고정 Blue→White→Red 스펙트럼 (원형 유지)"""
    norm_value = max(0.0, min(1.0, float(norm_value)))
    if norm_value <= 0.5:
        factor = norm_value * 2
        r = int(33 + (255 - 33) * factor)
        g = int(102 + (255 - 102) * factor)
        b = int(172 + (255 - 172) * factor)
    else:
        factor = (norm_value - 0.5) * 2
        r = int(255 - (255 - 178) * factor)
        g = int(255 - (255 - 24) * factor)
        b = int(255 - (255 - 43) * factor)
    r = max(0, min(255, r)); g = max(0, min(255, g)); b = max(0, min(255, b))
    return f"rgb({r}, {g}, {b})"

def validate_color_consistency(variety_value, global_min, global_max, expected_color, actual_color, theme, category_name):
    """색상 일치 검증 (원형 유지)"""
    if global_max == global_min:
        return True, "Uniform values - consistency guaranteed"
    norm_val = (float(variety_value) - float(global_min)) / (float(global_max) - float(global_min))
    norm_val = max(0, min(1, norm_val))
    recalc_color = generate_enhanced_gradient_color(norm_val, theme, category_name)
    is_consistent = (expected_color == actual_color == recalc_color)
    details = {
        'norm_val': norm_val,
        'expected': expected_color,
        'actual': actual_color,
        'recalculated': recalc_color,
        'consistent': is_consistent,
        'spectrum': 'Blue(#2166ac) → White(#ffffff) → Red(#b2182b)'
    }
    return is_consistent, details

def calculate_enhanced_trait_color_scale(phenotype_data, trait_name, category_name=None):
    """
    ✅ DB 의존 제거 버전
    - 글로벌 스케일: PHENO_DF 전체에서 trait 범위 계산(get_trait_full_range_data)
    - phenotype_data: 기존처럼 {'data': [...]} dict를 받을 수도 있고,
                      None이거나 형식이 다르면 PHENO_DF를 사용
    반환: { variety_id(str): {color, value, ...}, ... }
    """
    print(f"🎨 Enhanced color scale (PHENO_DF) 계산: trait={trait_name}, category={category_name}")

    # 입력 데이터프레임 결정
    if phenotype_data and isinstance(phenotype_data, dict) and 'data' in phenotype_data:
        df = pd.DataFrame(phenotype_data['data'])
    else:
        df = PHENO_DF.copy()

    if df.empty or trait_name not in df.columns or _ID_COL not in df.columns:
        print("   ❌ DataFrame 비었거나 필요한 컬럼(id/trait) 없음")
        return {}

    # 카테고리 자동 감지
    if not category_name:
        category_name = detect_trait_category(trait_name)

    trait_type = get_trait_type_from_db(trait_name)
    theme = get_category_color_theme(category_name)
    print(f"   🏷️ Trait 타입: {trait_type}")

    # NaN 제거
    clean_df = df[df[trait_name].notna()].copy()
    if clean_df.empty:
        print("   ❌ 유효 데이터 없음")
        return {}

    color_mapping = {}

    if trait_type == 'numeric':
        # Global 범위(전체 PHENO_DF) 확보
        full_range = get_trait_full_range_data(trait_name)
        if not full_range:
            print("   ❌ Global range 없음 → 색상 매핑 중단")
            return {}

        gmin, gmax = full_range['min'], full_range['max']
        gavg = full_range.get('avg', (gmin + gmax)/2)
        gcount = full_range.get('count', 0)

        vals = _coerce_numeric(clean_df[trait_name]).dropna()
        if vals.empty:
            return {}

        if gmin == gmax:
            mid_color = generate_enhanced_gradient_color(0.5, theme, category_name)
            for idx, v in vals.items():
                vid = str(clean_df.loc[idx, _ID_COL])
                color_mapping[vid] = {
                    'color': mid_color,
                    'value': float(v),
                    'normalized': 0.5,
                    'global_min': gmin, 'global_max': gmax,
                    'global_avg': gavg, 'global_count': gcount,
                    'trait_type': 'numeric'
                }
        else:
            span = (gmax - gmin)
            for idx, v in vals.items():
                vid = str(clean_df.loc[idx, _ID_COL])
                nv = max(0.0, min(1.0, (float(v) - gmin) / span))
                color = generate_enhanced_gradient_color(nv, theme, category_name)
                color_mapping[vid] = {
                    'color': color,
                    'value': float(v),
                    'normalized': nv,
                    'global_min': gmin, 'global_max': gmax,
                    'global_avg': gavg, 'global_count': gcount,
                    'trait_type': 'numeric'
                }

    elif trait_type in ('categorical', 'date'):
        unique_vals = clean_df[trait_name].astype(str).unique()
        palette = get_enhanced_seaborn_palette(category_name, len(unique_vals))
        # 최대 15개 팔레트만 지원하던 로직 그대로 유지
        for i, val in enumerate(unique_vals[:15]):
            color = palette[i % len(palette)]
            mask = clean_df[trait_name].astype(str) == val
            cnt = int(mask.sum())
            for idx in clean_df[mask].index:
                vid = str(clean_df.loc[idx, _ID_COL])
                color_mapping[vid] = {
                    'color': color,
                    'value': val,
                    'category_index': i,
                    'trait_type': trait_type,
                    'variety_count': cnt
                }

    print(f"   🌈 색상 매핑 완료: {len(color_mapping)}개 품종")
    return color_mapping


def pheno_all() -> pd.DataFrame:
    """전체 표현형 데이터프레임 반환 (사본)."""
    if PHENO_DF is None or PHENO_DF.empty:
        return pd.DataFrame()
    return PHENO_DF.copy()


def pheno_by_id(variety_id: Union[str, int]) -> pd.DataFrame:
    """단일 품종 id에 해당하는 행 반환 (사본)."""
    if PHENO_DF is None or PHENO_DF.empty or variety_id is None:
        return pd.DataFrame()
    return PHENO_DF[PHENO_DF[_ID_COL].astype(str) == str(variety_id)].copy()


def pheno_by_ids(variety_ids: Iterable[Union[str, int]]) -> pd.DataFrame:
    """여러 품종 id 행 일괄 반환 (사본)."""
    ids = [str(v) for v in (variety_ids or []) if v is not None]
    if not ids or PHENO_DF is None or PHENO_DF.empty:
        return pd.DataFrame()
    return PHENO_DF[PHENO_DF[_ID_COL].astype(str).isin(ids)].copy()


def pheno_rows_for_nodes(available_nodes: Optional[List[dict]]) -> pd.DataFrame:
    """
    기존 preload_phenotype_data()가 하던 '노드에서 it_number/id 추출'을
    그대로 따라, 해당 행들만 반환합니다. (preload 개념 자체는 불필요)
    """
    if not available_nodes:
        return pd.DataFrame()
    ids: List[str] = []
    for node in available_nodes:
        if isinstance(node, dict):
            itn = node.get("it_number") or node.get("id")
            if itn and str(itn) not in ["N/A", "None"]:
                ids.append(str(itn))
    if not ids:
        return pd.DataFrame()
    return pheno_by_ids(ids)


def build_pheno_content(trait: str, it_numbers: Optional[list] = None):
    df_all = pheno_all()
    if df_all.empty or trait not in df_all.columns:
        return html.Div("❌ No phenotype data or trait not found")

    # 전체 분포(연속형 가정) — 필요하면 타입 판정 함수로 분기
    s = pd.to_numeric(df_all[trait], errors="coerce").dropna()
    figs = []
    if not s.empty:
        figs.append(dcc.Graph(figure=px.histogram(s, x=trait, nbins=30, title=f"{trait} — distribution (all)")))

    # 선택 품종 하이라이트
    if it_numbers:
        sub = pheno_by_ids(it_numbers)
        if not sub.empty:
            vals = sub[[ "id", trait ]].copy()
            figs.append(html.Div([
                html.H5("Selected varieties"),
                html.Ul([html.Li(f"{r['id']}: {r[trait]}") for _, r in vals.iterrows()])
            ]))

    return html.Div(figs or [html.Div("No valid values")])


#========================================샘플추가모달






@callback(
    Output("additional-varieties-modal", "is_open"),
    [
        Input("open-additional-modal-btn", "n_clicks"),
        Input("addv-close-btn", "n_clicks"),
        Input("addv-apply-btn", "n_clicks"),
    ],
    State("additional-varieties-modal", "is_open"),
    prevent_initial_call=True,
)
def toggle_additional_varieties_modal(open_clicks, close_clicks, apply_clicks, is_open):
    ctx = dash.callback_context  # Dash >=2.0
    if not ctx.triggered:
        return is_open

    # Dash >=2.12 has ctx.triggered_id; fall back for older versions
    trigger_id = getattr(ctx, "triggered_id", None)
    if trigger_id is None:
        trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]

    if trigger_id == "open-additional-modal-btn":
        return True
    if trigger_id in ("addv-close-btn", "addv-apply-btn"):
        return False
    return is_open



#====수정완
# 1) 초기 로딩: 첫 20개 옵션만
@callback(
    Output("opt1-pedigree-variety-search", "options"),
    Input("opt1-pedigree-variety-search", "id"),
    prevent_initial_call=False
)
def init_pedigree_variety_search_options2(_):
    try:
        all_options = get_pedigree_variety_options()
        return all_options[:20]
    except Exception as e:
        print("opt1 init error:", e)
        return []

# 2) 검색 자동완성: 검색어 필터 + 현재 선택 유지(멀티)
@callback(
    [
        Output("opt1-pedigree-variety-search", "options", allow_duplicate=True),
        Output("opt1-pedigree-variety-search", "value", allow_duplicate=True),
    ],
    [
        Input("opt1-pedigree-variety-search", "search_value"),
        Input("opt1-pedigree-variety-search", "value"),
    ],
    prevent_initial_call=True,
)
def update_pedigree_variety_options2(search_value, current_value):
    """
    🔹 자동완성 입력이 사라지지 않도록 보완한 버전
    - search_value 입력 중에는 value를 업데이트하지 않음
    - value 변경(선택) 시에는 options을 그대로 유지
    """
    from dash import no_update

    cur_vals = _as_list(current_value)
    all_options = get_pedigree_variety_options()
    ctx = dash.callback_context
    if not ctx.triggered:
        return no_update, no_update

    trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]

    # 🧭 1️⃣ 검색 중일 때 (search_value 입력)
    if trigger_id == "opt1-pedigree-variety-search" and search_value is not None:
        needle = str(search_value).lower()
        filtered = [
            opt for opt in all_options
            if needle in str(opt.get("label", "")).lower()
            or needle in str(opt.get("value", "")).lower()
        ][:50]

        # 선택된 값은 항상 포함
        selected_values = set(cur_vals)
        for v in selected_values:
            if v not in [opt["value"] for opt in filtered]:
                found = next((opt for opt in all_options if opt["value"] == v), None)
                if found:
                    filtered.insert(0, found)

        return filtered, no_update

    # 🧭 2️⃣ 값 선택 이벤트일 때 (value 변경)
    if trigger_id == "opt1-pedigree-variety-search":
        return no_update, cur_vals

    return no_update, no_update


@callback(
    Output("addv-input-opt1", "data"),
    Input("opt1-pedigree-variety-search", "value"),
)
def opt1_value_to_store(vals):
    if not vals:
        return []
    return sorted([str(v) for v in (vals if isinstance(vals, list) else [vals])])

#==============수정완
@callback(
    Output('add-filter-btn_fix', 'style'),
    Input('filter-x_fix', 'value')
)
def toggle_filter_button_fix(selected_resources):
    return {'display': 'block'} if selected_resources else {'display': 'none'}

@callback(
    Output({'type': 'categorical-filter_fix', 'index': MATCH}, 'value'),
    Input({'type': 'categorical-filter_fix', 'index': MATCH}, 'value'),
    State('filter-x_fix', 'value'),
    State('filter-chain-store_fix', 'data'),
    State({'type': 'categorical-filter_fix', 'index': MATCH}, 'id'),
    prevent_initial_call=True,
)
def handle_categorical_selection_fix(selected_values, resource_type_value, filter_chain, component_id):
    """선택값이 비었고 categorical이면 해당 컬럼 전체값 자동 주입"""
    filter_index = component_id['index']
    try:
        current_filter = next(
            (f for f in (filter_chain or {}).get('filters', []) if f['id'] == filter_index),
            None
        )
        if (not selected_values) and current_filter and resource_type_value:
            if current_filter['type'] == 'categorical':
                df = get_db_data(resource_type_value)
                column = current_filter['column']
                all_values = sorted([str(v) for v in df[column].dropna().unique()])
                return all_values
    except Exception:
        pass
    return selected_values



# ----- (B) 필터 관리 (추가/제거) -----
@callback(
    Output('dynamic-filters_fix', 'children', allow_duplicate=True),
    Output('filter-chain-store_fix', 'data', allow_duplicate=True),
    Output('filter-modal_fix', 'is_open', allow_duplicate=True),
    Output('filter-type-selector_fix', 'options', allow_duplicate=True),

    Input('add-filter-btn_fix', 'n_clicks'),
    Input('confirm-filter-btn_fix', 'n_clicks'),
    Input('close-modal-btn_fix', 'n_clicks'),
    Input({'type': 'remove-filter_fix', 'index': ALL}, 'n_clicks'),
    Input('filter-x_fix', 'value'),

    State('filter-type-selector_fix', 'value'),
    State('dynamic-filters_fix', 'children'),
    State('filter-chain-store_fix', 'data'),
    prevent_initial_call=True
)
def manage_filters_fix(add_clicks, confirm_clicks, close_clicks, remove_clicks,
                       resource_type_value, selected_column, existing_filters, filter_chain):
    ctx_triggered = ctx.triggered_id
    if existing_filters is None:
        existing_filters = []
    if filter_chain is None:
        filter_chain = {'filters': [], 'results': []}

    # 리소스 타입 → 선택 가능한 컬럼 옵션
    try:
        if resource_type_value:
            df = get_db_data(resource_type_value)
            all_columns = get_available_columns(df)  # [{'label':..., 'value':...}, ...]
            used_columns = {f['column'] for f in filter_chain.get('filters', [])}
            available_options = [col for col in all_columns if col['value'] not in used_columns]
        else:
            available_options = []
    except Exception:
        available_options = []

    # 제거
    if isinstance(ctx_triggered, dict) and ctx_triggered.get('type') == 'remove-filter_fix':
        target_filter_id = ctx_triggered['index']

        # (외부) 필터 UI 매니저가 있다면 동기화
        try:
            opt2_filter_manager.remove_filter(f"filter_{target_filter_id}")
        except Exception:
            pass

        updated_filters, new_chain_filters = [], []
        # 기존 동적 필터 컴포넌트들의 index를 1..N으로 재배열
        for i, filter_component in enumerate(existing_filters, 1):
            if isinstance(filter_component, dict) and 'props' in filter_component:
                header_div = next(
                    (child for child in filter_component['props'].get('children', [])
                     if isinstance(child, dict) and 'props' in child and 'children' in child['props']
                     and isinstance(child['props']['children'], list) and len(child['props']['children']) > 1
                     and isinstance(child['props']['children'][1], dict)
                     and child['props']['children'][1].get('props', {}).get('id', {}).get('index') == target_filter_id),
                    None
                )
                if not header_div:
                    updated_component = update_filter_component_ids_fix(filter_component, i)
                    updated_filters.append(updated_component)
                    # 체인 갱신
                    cand = [f for f in filter_chain['filters'] if f['id'] != target_filter_id]
                    if cand:
                        nf = cand[0].copy()
                        nf['id'] = i
                        nf['sequence'] = i - 1
                        new_chain_filters.append(nf)

        filter_chain['filters'] = new_chain_filters
        return updated_filters, filter_chain, False, available_options

    # Add filter(+) → 모달 오픈
    if ctx_triggered == 'add-filter-btn_fix':
        return existing_filters, filter_chain, True, available_options

    # cancel
    if ctx_triggered == 'close-modal-btn_fix':
        return existing_filters, filter_chain, False, available_options

    # add(확인) → 실제 필터 추가
    if ctx_triggered == 'confirm-filter-btn_fix' and selected_column:
        df = get_db_data(resource_type_value)
        next_filter_id, filter_number = opt2_filter_manager.get_next_available_id()
        if next_filter_id is None:
            # 최대 개수 도달 등의 내부 정책 예외
            return existing_filters, filter_chain, False, available_options

        column_type = get_column_type_from_db(selected_column)
        # ⚠️ 새 컴포넌트의 내부 ID들도 모두 *_fix 로 생성되게 create_filter_component 수정 필요
        new_filter = create_filter_component_fix(filter_number, selected_column, column_type, df,
                                             id_suffix="_fix")
        initial_filter_values = get_initial_filter_values(df, column_type, selected_column)

        filter_chain['filters'].append({
            'id': filter_number,
            'column': selected_column,
            'type': column_type,
            'sequence': len(opt2_filter_manager.filter_sequence),
            'filter_values': initial_filter_values
        })
        return existing_filters + [new_filter], filter_chain, False, available_options

    # 기본 반환
    return existing_filters, filter_chain, False, available_options

def update_filter_component_ids_fix(component, new_index):
    """동적 필터 컴포넌트의 dict 구조에서 id.index만 교체 + *_fix 타입 유지"""
    if not isinstance(component, dict) or 'props' not in component:
        return component
    updated_component = component.copy()
    props = updated_component['props']
    if 'children' in props:
        props['children'] = [update_child_component_ids_fix(child, new_index) for child in props['children']]
    return updated_component

def update_child_component_ids_fix(child, new_index):
    if not isinstance(child, dict) or 'props' not in child:
        return child
    updated_child = child.copy()
    props = updated_child['props']
    if 'id' in props and isinstance(props['id'], dict):
        # 타입은 그대로 두되, index만 갱신
        t = props['id'].get('type')
        if t and t.endswith('_fix'):
            props['id'] = {'type': t, 'index': new_index}
        else:
            # 안전: 혹시 *_fix가 아니면 강제로 *_fix를 붙여줌
            props['id'] = {'type': f"{t}_fix", 'index': new_index}
    if 'children' in props:
        if isinstance(props['children'], list):
            props['children'] = [update_child_component_ids_fix(g, new_index) for g in props['children']]
        elif isinstance(props['children'], dict):
            props['children'] = update_child_component_ids_fix(props['children'], new_index)
    return updated_child


@callback(
    Output("opt2-pedi-checklist_fix", "data", allow_duplicate=True),
    Output("opt2-pedi-checklist_fix", "selected_row_ids", allow_duplicate=True),
    Output("opt2-last-table-data", "data", allow_duplicate=True),
    # 👉 당신이 원하신 트리거 유지(각 filter card의 변화)
    Input({'type': 'filter-result_fix', 'index': ALL}, "children"),
    # 👉 체인/리소스는 State 그대로 두되, 값이 없으면 일찍 반환
    Input('filter-chain-store_fix', 'data'),
    Input('filter-x_fix', 'value'),
    # 👉 직전 사용자 선택(전역 저장분)
    State("opt2-selected-ids", "data"),
    # 👉 이미 Current Selection에 담긴 개수(200 상한 계산)
    State("addv-input-multi", "data"),
    prevent_initial_call=True,
)
def opt2_apply_filters_fix(_filter_results_children, chain_data, resource_types,
                           persisted_ids, current_added):
    # 0) 방어
    if not resource_types:
        return [], [], []
    df = get_db_data(resource_types)
    if df.empty or not chain_data or not chain_data.get("filters"):
        return [], [], []

    # 1) 체인 적용
    filters = sorted(chain_data["filters"], key=lambda x: x["sequence"])
    filtered = df.copy()
    for f in filters:
        filtered = apply_stored_filter(filtered, f)

    ids = filtered["id"].tolist() if "id" in filtered.columns else []
    if not ids:
        return [], [], []

    # 2) 결과 합성 → 표준화
    r = build_results_table_data(ids)  # IT Number, Rice Variety Name, VCF Status, detection_status ...
    r = r.rename(columns={"IT Number": "id", "Rice Variety Name": "name"})
    if "VCF Status" in r.columns and "vcf_status" not in r.columns:
        r["vcf_status"] = r["VCF Status"]
    if "detection_status" not in r.columns:
        r["detection_status"] = "Not found"

    # 3) pedigree bool + 라벨(pedi/no pedi)
    r["has_pedigree"]   = r["detection_status"].apply(_pedigree_to_bool)
    r["pedigree_label"] = r["has_pedigree"].map(lambda b: "pedi" if b else "no pedi")

    out = r[["id", "name", "vcf_status", "pedigree_label", "has_pedigree"]].copy()
    out["id"] = out["id"].astype(str)

    # 4) 선택 복원: (현재 테이블 id) ∩ (persisted_ids)
    current_ids = set(out["id"])
    persisted   = set(map(str, persisted_ids or []))
    restored    = sorted(current_ids & persisted)

    # 5) 200 상한 적용(누적 기준): 이미 추가된 Current Selection 개수 반영
    cur_added = set(map(str, current_added or []))
    remaining = max(0, MAX_SEL - len(cur_added))
    if len(restored) > remaining:
        restored = restored[:remaining]

    # 6) 반환: data / selected_row_ids / last-table-data
    data_records = out.to_dict("records")
    return data_records, restored, data_records

@callback(
    Output("opt2-selected-ids", "data"),
    Input("opt2-pedi-checklist_fix", "selected_row_ids"),
)
def persist_selected_ids(selected_ids_now):
    return sorted(set(map(str, selected_ids_now or [])))

@callback(
    Output("opt2-pedi-checklist_fix", "selected_row_ids", allow_duplicate=True),
    Input("opt2-pedi-checklist_fix", "data"),
    Input("opt2-selected-ids", "data"),
    prevent_initial_call=True,
)
def sync_selection_to_table(rows, persisted):
    cur = {str(r.get("id")) for r in (rows or []) if r.get("id") is not None}
    per = set(map(str, (persisted or [])))
    return sorted(cur & per)

@callback(
    Output("opt2-selected-ids", "data", allow_duplicate=True),
    Input("addv-input-multi", "data"),
    State("opt2-last-table-data", "data"),
    prevent_initial_call=True,
)
def mirror_current_to_opt2_selected_ids(current_ids, last_rows):
    # 칩의 최종 상태를 테이블 선택에 미러링
    want = set(map(str, (current_ids or [])))
    df = pd.DataFrame(last_rows or [])
    if df.empty or "id" not in df.columns:
        return sorted(want)
    df["id"] = df["id"].astype(str)
    table_ids = set(df["id"].tolist())
    return sorted(want & table_ids)

# ----- (D) 체크 분배: 체크 → input-multiple, 나머지 → additional(phenotype) -----

@callback(
    Output('dynamic-filters_fix', 'style'),
    Input('filter-chain-store_fix', 'data')
)
def adjust_filter_container_height_fix(_):
    return {
        "maxHeight": "560px",
        "overflowY": "auto",
        "padding": "4px 0"
    }

def _safe_get(lst, ix, default=None):
    try:
        if lst is None: return default
        return lst[ix] if 0 <= ix < len(lst) else default
    except Exception:
        return default

@callback(
    Output('filter-chain-store_fix', 'data'),
    [
        Input({'type': 'filter-startdate-month_fix', 'index': ALL}, 'value'),
        Input({'type': 'filter-startdate-day_fix',   'index': ALL}, 'value'),
        Input({'type': 'filter-enddate-month_fix',   'index': ALL}, 'value'),
        Input({'type': 'filter-enddate-day_fix',     'index': ALL}, 'value'),
        Input({'type': 'numeric-range_fix',          'index': ALL}, 'value'),
        Input({'type': 'categorical-filter_fix',     'index': ALL}, 'value'),
    ],
    State('filter-chain-store_fix', 'data'),
    prevent_initial_call=True,
)
def update_filter_chain_fix(start_months, start_days, end_months, end_days,
                            range_values, categorical_values, current_chain):
    chain = current_chain or {'filters': [], 'results': []}
    trig = getattr(ctx, 'triggered_id', None)
    if not trig:
        return chain

    # 패턴 트리거만 처리 (dict 형태)
    if not isinstance(trig, dict):
        return chain

    ftype  = trig.get('type', '')
    findex = trig.get('index', None)
    if findex is None:
        return chain

    # 1-base index → 0-base로 변환
    i0 = int(findex) - 1

    # 대상 필터 조회
    target = next((f for f in chain.get('filters', []) if f.get('id') == findex), None)
    if not target:
        return chain

    # 날짜형
    if ftype.startswith('filter-startdate') or ftype.startswith('filter-enddate'):
        sm = _safe_get(start_months, i0)
        sd = _safe_get(start_days,   i0)
        em = _safe_get(end_months,   i0)
        ed = _safe_get(end_days,     i0)
        if all(v is not None for v in [sm, sd, em, ed]):
            try:
                start_date = int(sm) * 100 + int(sd)
                end_date   = int(em) * 100 + int(ed)
                target['filter_values'] = {
                    'date_range': [start_date, end_date],
                    'start_month': sm, 'start_day': sd,
                    'end_month': em,   'end_day': ed
                }
            except Exception:
                pass
        return chain

    # 수치형
    if ftype == 'numeric-range_fix':
        rv = _safe_get(range_values, i0)
        if rv:
            target['filter_values'] = {'range': rv}
        return chain

    # 범주형
    if ftype == 'categorical-filter_fix':
        cv = _safe_get(categorical_values, i0)
        if cv:
            target['filter_values'] = {'categories': cv}
        return chain

    return chain


@callback(
    Output("opt2-checklist-count", "children"),
    Input("opt2-pedi-checklist_fix", "data"),
)
def show_count(rows):
    n = len(rows or [])
    return f"Checklist (filtered results) — {n} varieties"

@callback(
    Output("addv-input-multi", "data", allow_duplicate=True),
    Output("opt2-notice", "children"),
    Input("opt2-selected-ids", "data"),
    State("opt2-last-table-data", "data"),
    State("addv-input-multi", "data"),
    prevent_initial_call=True,
)
def push_checked_to_current(selected_ids, last_rows, current):
    cur = set(map(str, current or []))
    sel = set(map(str, selected_ids or []))

    # ✅ 테이블에서 해제된 것은 칩에서도 제거
    kept = cur & sel

    remaining = MAX_SEL - len(kept)
    notice = None
    if remaining <= 0:
        if sel - kept:
            notice = dbc.Alert(
                f"Selection limit reached ({MAX_SEL}). No additional items were added.",
                color="warning", dismissable=True, duration=4000
            )
        return sorted(kept), notice

    # 새로 체크된 것만 추가
    new_candidates = sel - kept
    take = sorted(list(new_candidates))[:max(0, remaining)]
    merged = sorted(kept | set(take))

    skipped = len(new_candidates) - len(take)
    if skipped > 0:
        notice = dbc.Alert(
            f"Only {len(take)} item(s) added (limit {MAX_SEL}). {skipped} item(s) were skipped.",
            color="warning", dismissable=True, duration=4000
        )
    return merged, notice

@callback(
    Output("addv-additional-phenotype", "data"),
    Output("addv-additional-vcf", "data"),
    Input("addv-input-multi", "data"),
    State("opt2-last-table-data", "data"),
    prevent_initial_call=True,
)
def split_selection_meta(current_ids, last_rows):
    import re
    ids = set(map(str, current_ids or []))
    if not ids:
        return [], []

    df = pd.DataFrame(last_rows or [])
    if df.empty or "id" not in df.columns:
        return [], []

    df["id"] = df["id"].astype(str)

    # ✅ phenotype: pedigree_label == "pedi"
    pedi_series = df.get("pedigree_label", pd.Series([], dtype=str)).astype(str).str.lower()
    pheno_list = (
        df.loc[(df["id"].isin(ids)) & (pedi_series == "pedi"), "id"]
          .astype(str).drop_duplicates().sort_values().tolist()
    )

    # ✅ vcf: vcf_status/VCF Status 중 하나 사용 + RWG-#### 확실히 인식
    s = None
    if "vcf_status" in df.columns:
        s = df["vcf_status"].astype(str)
    elif "VCF Status" in df.columns:
        s = df["VCF Status"].astype(str)
    else:
        s = pd.Series([""] * len(df))

    s_clean = s.str.strip()
    has_rwg = s_clean.str.contains(r"(?i)\brwg-?\s*\d+", na=False)  # RWG-####
    not_empty = s_clean.ne("") & ~s_clean.str.fullmatch(r"(?i)\s*(no\s*vcf|none|null|nan)\s*")
    vcf_mask = has_rwg | not_empty

    vcf_list = (
        df.loc[(df["id"].isin(ids)) & vcf_mask, "id"]
          .astype(str).drop_duplicates().sort_values().tolist()
    )

    return pheno_list, vcf_list

@callback(
    Output("opt2-pedi-checklist_fix", "selected_row_ids", allow_duplicate=True),
    Input("opt2-pedi-checklist_fix", "selected_row_ids"),
    State("addv-input-multi", "data"),
    prevent_initial_call=True,
)
def hard_cap_selected_row_ids(selected_ids_now, current):
    cur = set(map(str, current or []))
    remaining = MAX_SEL - len(cur)
    s_now = [str(x) for x in (selected_ids_now or [])]
    if remaining < len(s_now):
        # 남은 슬롯만 허용
        return s_now[:max(0, remaining)]
    return dash.no_update

@callback(
    Output("opt2-quick-panel", "children"),
    Input("opt2-selected-ids", "data"),
    State("opt2-last-table-data", "data"),
)
def render_chips(current_ids, last_rows):
    ids = [str(x) for x in (current_ids or [])]
    if not ids:
        return [html.Span("No items selected yet.", style={"color": "#adb5bd"})]

    import pandas as pd
    df = pd.DataFrame(last_rows or [])
    if df.empty:
        # 데이터프레임이 비어도 칩은 보여줍니다 (IT만)
        return [
            html.Span(
                [
                    html.Span(i),
                    html.Button("×", id={"type": "opt2-chip-remove", "id": i}, n_clicks=0, style=chip_x_style),
                ],
                style=chip_style,
            )
            for i in ids[:15]
        ]

    # id는 항상 문자열
    df["id"] = df.get("id", pd.Series(dtype=str)).astype(str)

    # name은 원시값 유지(astype로 "None" 같은 문자열화 금지)
    nm_raw = df["name"] if "name" in df.columns else pd.Series([None] * len(df), index=df.index)
    name_map = dict(zip(df["id"], nm_raw))

    # ✅ None/"", "none"/"nan"/"null" → 이름 없음으로 간주 → IT만 표시
    def _label(i: str) -> str:
        n = name_map.get(i)
        if n is None:
            return i
        n = str(n).strip()
        if n == "" or n.lower() in {"none", "nan", "null"}:
            return i
        return f"{n} ({i})"

    cap = 15
    chips = []
    for vid in ids[:cap]:
        chips.append(
            html.Span(
                [
                    html.Span(_label(vid)),
                    html.Button("×", id={"type": "opt2-chip-remove", "id": vid}, n_clicks=0, style=chip_x_style),
                ],
                style=chip_style,
            )
        )
    if len(ids) > cap:
        chips.append(html.Span(f"… +{len(ids)-cap} more", style={"color": "#868e96"}))
    return chips

@callback(
    Output("addv-input-multi", "data", allow_duplicate=True),
    Input({"type":"opt2-chip-remove","id": ALL}, "n_clicks"),
    State("addv-input-multi", "data"),
    State({"type":"opt2-chip-remove","id": ALL}, "id"),
    prevent_initial_call=True,
)
def remove_by_chip(clicks, current, ids_json):
    if not clicks or not any(clicks):
        return dash.no_update
    rm = {str(ids_json[i]["id"]) for i, c in enumerate(clicks) if c}
    return sorted([x for x in (current or []) if str(x) not in rm])

@callback(
    Output("opt2-selected-ids", "data", allow_duplicate=True),
    Input({"type":"opt2-chip-remove","id": ALL}, "n_clicks"),
    State("opt2-selected-ids", "data"),
    State({"type":"opt2-chip-remove","id": ALL}, "id"),
    prevent_initial_call=True,
)
def remove_by_chip_opt2(clicks, current, ids_json):
    if not clicks or not any(clicks):
        return dash.no_update
    rm = {str(ids_json[i]["id"]) for i, c in enumerate(clicks) if c}
    return sorted([x for x in (current or []) if str(x) not in rm])


@callback(
    Output("opt2-quick-input-fix", "options"),
    Input("opt2-last-table-data", "data"),
)
def fill_quick_options(rows):
    df = pd.DataFrame(rows or [])
    if df.empty:
        return []
    return [{"label": f'{r.get("name") or r["id"]} ({r["id"]})', "value": str(r["id"])} for _, r in df.iterrows()]

@callback(
    Output("addv-input-multi", "data", allow_duplicate=True),
    Input("opt2-quick-input-fix", "value"),
    State("addv-input-multi", "data"),
    prevent_initial_call=True,
)
def quick_add(vals, cur):
    cur = set(cur or [])
    for v in (vals or []):
        cur.add(str(v))
    return sorted(cur)

@callback(
    Output("opt2-quick-input-fix", "value", allow_duplicate=True),
    Input("addv-input-multi", "data"),
    prevent_initial_call=True,
)
def mirror_current_to_quick(current_ids):
    return sorted([str(x) for x in (current_ids or [])])

@callback(
    Output("opt2-selected-ids", "data", allow_duplicate=True),
    Input("opt2-quick-input-fix", "value"),
    prevent_initial_call=True,
)
def quick_to_opt2_selected_ids(vals):
    return sorted([str(v) for v in (vals or [])])

@callback(
    Output("addv-input-opt2", "data"),
    Input("addv-input-multi", "data"),
)
def mirror_opt2_to_store(v):
    return sorted([str(x) for x in (v or [])])

#=========================
def _pill_style(active=False):
    base = {"padding":"4px 8px","borderRadius":"999px","border":"1px solid #d0d7de",
            "cursor":"pointer","userSelect":"none"}
    if active:
        base.update({"background":"#e8f0fe","border":"1px solid #90caf9","fontWeight":"600"})
    return base

@callback(
    Output('ecotype-filter-container2', 'children'),
    Output('variety-group-filter-container2', 'children'),
    Input('vcf-info-global', 'data'),   # 있으면 사용, 없으면 VCF_INFO_DF
    Input('opt3-active-filter', 'data'),
)
def render_filter_pills(vcf_info_data, active):
    vcf_df = pd.DataFrame(vcf_info_data) if vcf_info_data else VCF_INFO_DF.copy()
    if vcf_df.empty:
        return [], []

    ecos = [e for e in vcf_df.get('Ecotype', pd.Series(dtype=str)).dropna().unique()]
    main = ['Temperate japonica','Indica','Tropical japonica']
    other = sorted([e for e in ecos if e not in main])

    eco_children=[]
    for ec in main:
        if ec in ecos:
            flag = (active and active.get('type')=='ecotype' and active.get('value')==ec)
            eco_children.append(html.Span(ec, id={'type':'opt3-eco','value':ec},
                                          className='ecotype-span-hover', style=_pill_style(flag)))
    if other:
        flag = (active and active.get('type')=='ecotype' and active.get('value')=='Other')
        eco_children.append(html.Span("Other", id={'type':'opt3-eco','value':'Other'},
                                      className='ecotype-span-hover', style=_pill_style(flag)))

    vgs = sorted([v for v in vcf_df.get('Variety Group', pd.Series(dtype=str)).dropna().unique()])
    vg_children=[]
    for vg in vgs:
        flag = (active and active.get('type')=='variety_group' and active.get('value')==vg)
        vg_children.append(html.Span(vg, id={'type':'opt3-vg','value':vg},
                                     className='variety-group-span-hover', style=_pill_style(flag)))
    return eco_children, vg_children

@callback(
    Output("opt3-active-filter", "data"),
    Input({'type':'opt3-eco','value':ALL}, 'n_clicks'),
    Input({'type':'opt3-vg','value':ALL}, 'n_clicks'),
    prevent_initial_call=True,
)
def click_opt3_pills(ec_clicks, vg_clicks):
    trig = getattr(dash.callback_context, "triggered_id", None)
    if isinstance(trig, dict):
        if trig.get("type")=="opt3-eco":
            return {"type":"ecotype","value":trig.get("value")}
        if trig.get("type")=="opt3-vg":
            return {"type":"variety_group","value":trig.get("value")}
    return dash.no_update

# 1️⃣ 초기 로딩: 첫 20개 옵션만
@callback(
    Output("opt3-vcf-sample-dropdown", "options"),
    Input("opt3-vcf-sample-dropdown", "id"),
    prevent_initial_call=False
)
def init_vcf_sample_options2(_):
    try:
        df = VCF_INFO_DF.copy()
        if df.empty:
            return []
        req = ["Variety Name in English", "Entry_No."]
        if any(c not in df.columns for c in req):
            return []

        opts = []
        for _, r in df.iterrows():
            name = str(r["Variety Name in English"]).strip()
            entry = str(r["Entry_No."]).strip()
            group = str(r.get("Variety Group", "")).strip()
            label = f"{name} ({entry})" + (f" · {group}" if group else "")
            opts.append({"label": label, "value": entry})

        return opts[:20]
    except Exception as e:
        print("opt3 init error:", e)
        return []

@callback(
    [
        Output("opt3-vcf-sample-dropdown", "options", allow_duplicate=True),
        Output("opt3-vcf-sample-dropdown", "value", allow_duplicate=True),
    ],
    [
        Input("opt3-vcf-sample-dropdown", "search_value"),
        Input("opt3-vcf-sample-dropdown", "value"),
    ],
    prevent_initial_call=True,
)
def update_vcf_sample_options(search_value, current_value):
    """
    🔹 opt3 자동완성 기능
    - 여러 항목 선택 가능
    - checklist와의 sync 방해하지 않음
    """
    from dash import no_update
    import dash

    cur_vals = _as_list(current_value)
    ctx = dash.callback_context
    if not ctx.triggered:
        return no_update, no_update

    trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]

    df = VCF_INFO_DF.copy()
    if df.empty:
        return [], cur_vals

    # 전체 옵션 구성
    opts_all = []
    for _, r in df.iterrows():
        name = str(r["Variety Name in English"]).strip()
        entry = str(r["Entry_No."]).strip()
        group = str(r.get("Variety Group", "")).strip()
        label = f"{name} ({entry})" + (f" · {group}" if group else "")
        opts_all.append({"label": label, "value": entry})

    # 1️⃣ search_value 입력 중
    if trigger_id == "opt3-vcf-sample-dropdown" and search_value is not None:
        needle = str(search_value).lower()
        filtered = [
            opt for opt in opts_all
            if needle in str(opt.get("label", "")).lower()
            or needle in str(opt.get("value", "")).lower()
        ][:50]

        # 선택값 유지
        existing_vals = set(cur_vals)
        for v in existing_vals:
            if v not in [opt["value"] for opt in filtered]:
                found = next((opt for opt in opts_all if opt["value"] == v), None)
                if found:
                    filtered.insert(0, found)

        return filtered, no_update

    # 2️⃣ 값 선택 이벤트 시
    if trigger_id == "opt3-vcf-sample-dropdown":
        return no_update, cur_vals

    return no_update, no_update

@callback(
    Output("opt3-vcf-checklist", "data"),
    Input("opt3-active-filter", "data"),
)
def build_checklist_by_pills(active_filter):
    df = VCF_INFO_DF.copy()
    if df.empty or not active_filter:
        return []
    t, v = active_filter.get("type"), active_filter.get("value")
    if t == "ecotype":
        if v == "Other":
            main = {'Temperate japonica','Indica','Tropical japonica'}
            df = df[~df["Ecotype"].isin(main)]
        else:
            df = df[df["Ecotype"]==v]
    elif t == "variety_group":
        df = df[df["Variety Group"]==v]

    # variety_id 매핑(있으면 채움) — load_vcf_mapping 사용
    try:
        vmap = load_vcf_mapping()  # columns: variety_id, vcf_tag(=Entry_No.)
        vmap.rename(columns={"vcf_tag":"Entry_No."}, inplace=True)  # 통일
        df = df.merge(vmap, on="Entry_No.", how="left")
    except Exception:
        df["variety_id"] = ""

    cols = ["Entry_No.","Variety Name in English","Ecotype","Variety Group","variety_id"]
    cols = [c for c in cols if c in df.columns]
    return df[cols].to_dict("records")

def _variety_has_pedigree(it_id: str) -> bool:
    if not it_id: return False
    row = VARIETY_DF[VARIETY_DF["id"]==it_id]
    return (not row.empty) and bool(row.iloc[0].get("has_pedigree"))

@callback(
    Output("opt3-vcf-sample-dropdown", "value", allow_duplicate=True),
    Output("opt3-vcf-checklist", "selected_rows", allow_duplicate=True),
    Input("opt3-vcf-checklist", "selected_rows"),
    State("opt3-vcf-checklist", "data"),
    State("opt3-vcf-sample-dropdown", "value"),
    prevent_initial_call=True,
)
def checklist_to_dropdown(selected_rows, checklist_data, current_values):
    import pandas as pd
    df = pd.DataFrame(checklist_data or [])
    vals = set(map(str, (current_values or [])))
    if df.empty:
        return sorted(vals), dash.no_update

    # 현재 화면(필터/페이지)에서 보이는 모든 Entry_No.
    visible = set(df.get("Entry_No.", pd.Series(dtype=str)).astype(str).tolist())

    # 체크된 항목
    picked = set()
    if selected_rows:
        picked = set(df.iloc[selected_rows]["Entry_No."].astype(str).tolist())

    # 보이는 항목 중 해제된 것들은 드롭다운에서 제거,
    # 체크된 것은 추가 (보이지 않는 항목은 유지)
    new_vals = (vals - (visible - picked)) | picked

    return sorted(new_vals), dash.no_update

@callback(
    Output("opt3-vcf-checklist", "selected_rows", allow_duplicate=True),
    Input("opt3-vcf-sample-dropdown", "value"),
    Input("opt3-vcf-checklist", "data"),
    prevent_initial_call=True,
)
def sync_checklist_selection(drop_values, checklist_data):
    import pandas as pd
    df = pd.DataFrame(checklist_data or [])
    if df.empty:
        return []

    want = set(map(str, (drop_values or [])))
    if not want:
        return []

    # 현재 표시 중인 rows 중, 드롭다운 값과 일치하는 행의 인덱스 선택
    sr = df.index[df.get("Entry_No.", pd.Series(dtype=str)).astype(str).isin(want)].tolist()
    return sr



@callback(
    Output("opt3-preview-table", "data"),
    Output("opt3-selected-preview-rows", "data"),
    Output("addv-input-opt3", "data"),
    Input("opt3-vcf-sample-dropdown", "value"),
)
def update_opt3_preview_from_dropdown(selected_ids):
    import pandas as pd
    df = VCF_INFO_DF.copy()
    if df.empty or not selected_ids:
        return [], [], []
    want = [str(x) for x in (selected_ids if isinstance(selected_ids, list) else [selected_ids])]
    df = df[df["Entry_No."].astype(str).isin(want)]

    # ✅ 프리뷰/Store는 이 3컬럼만 확정 사용
    cols3 = ["Entry_No.","Variety Name in English","Variety Group","Ecotype","Origin","NCBI Biosample ID"]

    have = [c for c in cols3 if c in df.columns]
    data = df[have].to_dict("records")
    print(data)

    entries = sorted(df["Entry_No."].astype(str).unique().tolist())
    return data, data, entries



@callback(
    [
#Output("addv-apply-result", "children"),
        Output("additional-varieties-modal", "is_open", allow_duplicate=True),
        Output("addv-apply-groups-pedi", "data"),
        Output("addv-apply-groups-nopedi", "data"),
        Output("addv-pedi-message", "children"),  # ✅ 추가: 안내 메시지 출력용
    ],
    Input("addv-apply-btn", "n_clicks"),
    State("addv-meta-cache", "data"),
    prevent_initial_call=True,
)
def on_apply(n, meta_cache):
    import pandas as pd, re
    from dash import html
    from dash.exceptions import PreventUpdate

    df_all = pd.DataFrame(meta_cache or [])
    if df_all.empty:
        pedi_message = html.Div([
            html.Span("👉 Use "),
        html.B('"Add Additional Varieties"'),
        html.Span(" to include more varieties and start analysis."),
        ], style={
            'background': '#fff3cd',
            'padding': '10px',
            'borderRadius': '6px',
            'fontSize': '0.9rem'
        })

        group_pedi = []  # ✅ 빈 리스트
        group_nopedi = {"records": []}  # ✅ 빈 dict 구조 유지

        # 필요한 output 구조에 맞게 False 또는 no_update 처리
        return False, group_pedi, group_nopedi, pedi_message
    # ====================================================
    # 🟩 1. PEDI 그룹: 계보 탐색 가능 품종들
    # ====================================================
    pedi_df = df_all[df_all["pedigree_label"].str.lower() == "pedi"].copy()
    group_pedi = sorted(pedi_df["name"].dropna().astype(str).unique().tolist())

    # ====================================================
    # 🟥 2. NOPEDI 그룹: 계보 탐색 불가 품종들
    # ====================================================
    nopedi_df = df_all[df_all["pedigree_label"].str.lower() != "pedi"].copy()
    nopedi_records = []

    for _, row in nopedi_df.iterrows():
        name = str(row.get("name", "")).strip()
        it_number = str(row.get("id", "")).strip()
        vcf_status = str(row.get("vcf_status", "")).strip()
        pedi_label = str(row.get("pedigree_label", "no pedi")).strip().lower() or "no pedi"

        has_vcf = bool(vcf_status and not re.match(r"(?i)^\s*(no\s*vcf|none|null|nan)\s*$", vcf_status))
        has_pheno = bool(it_number)
        vcf_id = vcf_status if has_vcf else ""

        nopedi_records.append({
            "id": name,
            "name": name,
            "it_number": it_number,
            "has_vcf": has_vcf,
            "vcf_id": vcf_id,
            "has_pheno": has_pheno,
            "pedigree_label": pedi_label,
        })

    group_nopedi = {"records": nopedi_records}
    #print(group_nopedi)

    # ====================================================
    # 🟨 3. 결과 메시지 처리
    # ====================================================
    total = len(df_all)
    n_pedi = len(group_pedi)
    n_nopedi = len(nopedi_records)

    txt = f"[Apply] 🟩 pedi: {n_pedi} | 🟥 nopedi: {n_nopedi} → total: {total}"

    # --- 안내 메시지 구성
    if n_pedi == 0:
        pedi_message = html.Div([
            html.Span("❌ No pedigree-available varieties were found among your selection.", 
                      style={'color': '#dc3545', 'fontWeight': '600'}),
            html.Br(),
            html.Span("You can still visualize the non-pedigree varieties, or click "),
            html.B('"Add Additional Varieties"'),
            html.Span(" to explore varieties with full pedigree information."),
        ], style={'background': '#fff3cd', 'padding': '8px', 'borderRadius': '6px', 'fontSize': '0.9rem'})
    else:
        pedi_message = html.Div([
            html.Span("✅ Pedigree visualization is available for the selected varieties.", 
                      style={'color': '#28a745', 'fontWeight': '600'}),
            html.Br(),
            html.Span(f"Found {n_pedi} variety{'ies' if n_pedi > 1 else ''} with pedigree data ready for analysis.")
        ], style={'background': '#e8f5e9', 'padding': '8px', 'borderRadius': '6px', 'fontSize': '0.9rem'})

    #print("🚩 group_pedi:", len(group_pedi))
    #print("🚩 group_nopedi (sample):", nopedi_records[:3])

    return False, group_pedi, group_nopedi, pedi_message




  #======== 마지막정리




# ✅ 추가: opt3 드롭다운 Entry_No. → IT id로 매핑해 저장


# ✅ 유틸(파일 상단 또는 utils 모듈)
def _meta_for_ids(ids):
    import pandas as pd, numpy as np, re
    ids = sorted(set(str(x) for x in (ids or [])))
    if not ids:
        return pd.DataFrame(columns=["id","name","pedigree_label","vcf_status"])

    r = build_results_table_data(ids)  # IT Number, Rice Variety Name, VCF Status, detection_status ...
    if r is None or r.empty:
        return pd.DataFrame(columns=["id","name","pedigree_label","vcf_status"])

    r = r.rename(columns={"IT Number":"id","Rice Variety Name":"name"})
    if "VCF Status" in r.columns and "vcf_status" not in r.columns:
        r["vcf_status"] = r["VCF Status"]
    r["id"] = r["id"].astype(str)

    # --- pedi 신호들: 전부 r.index 길이로 맞춤 ---
    # A) detection_status -> found == pedi
    if "detection_status" in r.columns:
        ds_pedi = r["detection_status"].astype(str).str.lower().eq("found").fillna(False)
    else:
        ds_pedi = pd.Series(False, index=r.index)

    # B) 기존 pedigree_label == 'pedi'
    if "pedigree_label" in r.columns:
        pl_pedi = r["pedigree_label"].astype(str).str.lower().eq("pedi").fillna(False)
    else:
        pl_pedi = pd.Series(False, index=r.index)

    # C) VARIETY_DF.has_pedigree 를 map으로 (merge 하지 않음)
    try:
        vd = VARIETY_DF[["id","has_pedigree"]].copy()
        vd["id"] = vd["id"].astype(str)
        vd_map = vd.set_index("id")["has_pedigree"].astype(bool)
        vd_pedi = r["id"].map(vd_map).fillna(False)
    except Exception:
        vd_pedi = pd.Series(False, index=r.index)

    # 최종 판정: 하나라도 True면 pedi
    any_pedi = (ds_pedi | pl_pedi | vd_pedi).astype(bool)
    r["pedigree_label"] = np.where(any_pedi, "pedi", "no pedi")

    # vcf_status 정규화
    if "vcf_status" in r.columns:
        r["vcf_status"] = r["vcf_status"].fillna("NaN").astype(str)
        r.loc[
            r["vcf_status"].str.fullmatch(r"(?i)\s*(nan|none|null)\s*"),
            "vcf_status"
        ] = ""

    return r[["id","name","pedigree_label","vcf_status"]].drop_duplicates()

@callback(
    Output("addv-union", "data"),
    Output("addv-meta-cache", "data"),
    Input("addv-input-opt1", "data"),
    Input("addv-input-opt2", "data"),
    Input("addv-input-opt3", "data"),
    Input("opt3-selected-preview-rows", "data"),
)
def build_union_and_meta(opt1_ids, opt2_ids, opt3_ids, opt3_preview_rows):
    import pandas as pd, re, unicodedata

    # ---- 공통 헬퍼 ----
    def _safe_str(s): return str(s) if s is not None and pd.notna(s) else ""
    def _looks_like_it(s): return bool(re.match(r"^IT\d+$", str(s).strip(), re.I))
    def _looks_like_entry(s): return bool(re.match(r"^[A-Za-z]{2,5}-?\d+$", str(s).strip()))
    def _normalize_text_no_space(s): return re.sub(r"\s+", "", _safe_str(s).lower())
    def _compact_line_name(s): return re.sub(r"\s+", "", re.sub(r"(\d+)\s+", r"\1", _safe_str(s).lower()))
    def _match_in_thal3(name, th_names):
        if not name: return False
        nrm = unicodedata.normalize("NFC", name.strip().lower().replace("-", "").replace(" ", ""))
        for t in th_names:
            tn = unicodedata.normalize("NFC", t.strip().lower().replace("-", "").replace(" ", ""))
            if nrm == tn: return True
        return False
    print(opt3_preview_rows)
    # ============================================================
    # 1️⃣ 각 입력 세트별로 구분
    # ============================================================
    opt1_ids = [str(x) for x in (opt1_ids or [])]
    opt2_ids = [str(x) for x in (opt2_ids or [])]
    opt3_ids = [str(x) for x in (opt3_ids or [])]

    # 전체 union (나중에 출력용)
    union_ids = sorted(set(opt1_ids + opt2_ids + opt3_ids))

    # ============================================================
    # 2️⃣ opt2 : IT 기반
    # ============================================================
    df_it = _meta_for_ids([i for i in opt2_ids if _looks_like_it(i)])
    df_it = df_it[["id","name","pedigree_label","vcf_status"]] if not df_it.empty else pd.DataFrame(columns=["id","name","pedigree_label","vcf_status"])
    df_it["source"] = "opt2"

    # ============================================================
    # 3️⃣ opt3 : Entry_No. 기반
    # ============================================================
    
    df_en = pd.DataFrame(columns=["id", "name", "pedigree_label", "vcf_status", "source"])
    try:
        prev = pd.DataFrame(opt3_preview_rows or [])
        if not prev.empty:
            prev = prev.rename(columns={
                "Entry_No.": "Entry_No.",  # 그대로 유지 (안전)
                "Variety Name in English": "Variety (EN)",
                "variety_id": "variety_id",
            })

            for col in ["Entry_No.", "Variety (EN)", "variety_id"]:
                if col in prev.columns:
                    prev[col] = prev[col].astype(str)
                else:
                    prev[col] = ""

        # ✅ 1단계: vmap으로 Entry_No. → Variety_ID 매핑
        vmap = pd.DataFrame()
        try:
            vmap = load_vcf_mapping()  # columns: variety_id, vcf_tag(=Entry_No.)
            vmap.rename(columns={"vcf_tag": "Entry_No."}, inplace=True)
        except Exception as e:
            print(f"⚠️ load_vcf_mapping 실패: {e}")

        # ✅ vmap 기반 variety_id 보강
        if not vmap.empty:
            prev = prev.merge(vmap, on="Entry_No.", how="outer", suffixes=("", "_vmap"))
            prev["variety_id"] = prev["variety_id"].fillna(prev.get("variety_id_vmap", "")).fillna("")

        # ✅ 2단계: 계보 탐지용 reference 준비
        vd_has = {}
        try:
            vd_has = {str(k): bool(v) for k, v in zip(VARIETY_DF["id"], VARIETY_DF["has_pedigree"])}
        except Exception:
            pass

        th_names = set()
        try:
            th = load_thal3()
            cand_col = next((c for c in ["Variety (EN)", "cultivar_name"] if c in th.columns), None)
            if cand_col:
                th_names = set(th[cand_col].astype(str).str.strip())
        except Exception:
            pass

        # ✅ 3단계: opt3_ids 순회하며 Entry_No. 기반 최종 병합
        rows = []
        for eid in opt3_ids:  # eid == RWG-코드
            r = prev.loc[prev["Entry_No."] == eid]
            vname, vid = "", ""

            if not r.empty:
                rr = r.iloc[0]
                vname = str(rr.get("Variety (EN)", "")).strip()
                vid = str(rr.get("variety_id", "")).strip()

            # ✅ vmap fallback — prev에 없거나 vid 비어있을 때
            if not vid and not vmap.empty:
                match = vmap.loc[vmap["Entry_No."] == eid]
                if not match.empty:
                    vid = str(match.iloc[0].get("variety_id", "")).strip()

            # ✅ 계보 여부 판단
            pedi = "pedi" if (vid and vd_has.get(vid)) or _match_in_thal3(vname, th_names) else "no pedi"

            # ✅ 최종 row 추가
            rows.append({
                "id": vid,                           # ✅ Variety ID (IT number)
                "name": vname if vname else eid,     # ✅ 품종명 없으면 Entry_No.
                "pedigree_label": pedi,
                "vcf_status": eid,                   # ✅ Entry_No. (RWG code)
                "source": "opt3"
            })

        df_en = pd.DataFrame(rows)
    except Exception as e:
        print(f"⚠️ opt3 처리 오류: {e}")
    except Exception as e:
        print(f"⚠️ opt3 처리 오류: {e}")

    # ============================================================
    # 4️⃣ opt1 : 품종명 기반 (pedigree)
    # ============================================================
    df_nm = pd.DataFrame(columns=["id","name","pedigree_label","vcf_status"])
    try:
        rt = load_resource_types()
        vm = load_vcf_mapping()
        th = load_thal3()

        rt["variety_key"] = rt["variety_name_en"].map(_normalize_text_no_space)
        rt["line_key"] = rt["line_name_en"].map(_compact_line_name)
        vm_dict = dict(zip(vm["variety_id"].astype(str), vm["vcf_tag"].astype(str)))

        th["cultivar_key"] = th["cultivar_name"].map(_normalize_text_no_space)
        th["parent_key"] = th["parent"].map(_normalize_text_no_space)

        rows = []
        for nm in opt1_ids:
            nkey = _normalize_text_no_space(nm)
            it_id = ""
            vcf_tag = ""

            hit = rt.loc[(rt["variety_key"] == nkey) | (rt["line_key"] == nkey)]
            if not hit.empty:
                it_id = str(hit.iloc[0]["id"])
                vcf_tag = vm_dict.get(it_id, "")
            else:
                th_hit = th.loc[(th["cultivar_key"] == nkey) | (th["parent_key"] == nkey)]
                if not th_hit.empty:
                    cand = th_hit.iloc[0]["cultivar_name"]
                    cand_key = _normalize_text_no_space(cand)
                    hit2 = rt.loc[(rt["variety_key"] == cand_key) | (rt["line_key"] == cand_key)]
                    if not hit2.empty:
                        it_id = str(hit2.iloc[0]["id"])
                        vcf_tag = vm_dict.get(it_id, "")

            rows.append({"id": it_id, "name": nm, "pedigree_label": "pedi", "vcf_status": vcf_tag,'source':'opt1'})
        df_nm = pd.DataFrame(rows)
    except Exception as e:
        print(f"⚠️ opt1 처리 오류: {e}")

    # ============================================================
    # 5️⃣ 통합 및 반환
    # ============================================================
    df_all = pd.concat([df_nm, df_it, df_en], ignore_index=True).drop_duplicates(subset=["id","name"], keep="first")
    if not df_all.empty:
        df_all["vcf_status"] = df_all["vcf_status"].fillna("").astype(str)
        df_all.loc[df_all["vcf_status"].str.fullmatch(r"(?i)\s*(nan|none|null)\s*"), "vcf_status"] = ""

    return union_ids, df_all.to_dict("records")

@callback(
    Output("addv-input-multi-preview", "children"),
    Output("addv-additional-phenotype-preview", "children"),
    Output("addv-additional-vcf-preview", "children"),
    Output("addv-additional-pedigree-count", "children"),
    Output("addv-additional-phenotype-count", "children"),
    Output("addv-additional-vcf-count", "children"),
    Input("addv-union", "data"),
    State("addv-meta-cache", "data"),
)
def render_summary(union_ids, meta_rows):
    import pandas as pd, numpy as np
    from dash import html

    df = pd.DataFrame(meta_rows or [])
    if df.empty:
        return html.Div("No data available"), [], [], "Pedigree: 0", "Phenotype: 0", "VCF: 0"

    # ✅ 기본 전처리
    for col in ["id", "name", "pedigree_label", "vcf_status", "source"]:
        if col not in df.columns:
            df[col] = ""
    df = df.astype(str)

    # ✅ 색상 규칙 (opt1=초록, opt2=파랑, opt3=보라)
    COLOR_MAP = {"opt1": "#27ae60", "opt2": "#2980b9", "opt3": "#8e44ad"}

    # ✅ 칩 생성 함수 (source별 prefix / 색상 / 테두리)
    def make_chip(label, source, pedi_label):
        color = COLOR_MAP.get(source, "#7f8c8d")
        border_color = "#2ecc71" if pedi_label.lower() == "pedi" else "#e74c3c"
        prefix = {"opt1": "①", "opt2": "②", "opt3": "③"}.get(source, "")
        return html.Span(
            f"{prefix} {label}",
            style={
                "display": "inline-block",
                "padding": "3px 6px",
                "margin": "3px",
                "border": f"1px solid {border_color}",
                "borderRadius": "10px",
                "fontSize": "12px",
                "color": color,
            },
        )

    # ✅ source별 분류
    opt1_df = df[df["source"] == "opt1"].copy()
    opt2_df = df[df["source"] == "opt2"].copy()
    opt3_df = df[df["source"] == "opt3"].copy()

    # ✅ label 규칙 (id=IT번호, name=표시용 품종명)
    def fmt_opt1(row):
        return f"{row['name']} ({row['id']})" if row["id"] else row["name"]

    def fmt_opt2(row):
        name = row.get("name")
        vid = row.get("id")
        return f"{name} ({vid})" if name and vid else (name or str(vid) or None)

    def fmt_opt3(row):
        vc = row.get("vcf_status", "")
        return f"{row['name']} ({vc})" if vc else row["name"]

    opt1_df["label"] = opt1_df.apply(fmt_opt1, axis=1)
    opt2_df["label"] = opt2_df.apply(fmt_opt2, axis=1)
    opt3_df["label"] = opt3_df.apply(fmt_opt3, axis=1)

    # ✅ 칩 생성 (source별 순서 유지)
    chips = []
    for src_df in [opt1_df, opt2_df, opt3_df]:
        for _, row in src_df.iterrows():
            chips.append(make_chip(row["label"], row["source"], row.get("pedigree_label", "")))

    if len(chips) > 60:
        chips = chips[:60] + [html.Span(f" … +{len(df)-60} more", style={"color": "#95a5a6"})]

    # ✅ Count 계산
    n_pedi = (df["pedigree_label"].str.lower() == "pedi").sum()
    n_pheno = df["id"].replace("", np.nan).notna().sum()
    n_vcf = df["vcf_status"].replace("", np.nan).notna().sum()

    # ✅ Preview Lists
    pheno_list = [html.Div(fmt_opt2(r)) for _, r in opt2_df.head(200).iterrows()]
    vcf_list = [html.Div(fmt_opt3(r)) for _, r in opt3_df.head(200).iterrows()]

    return (
        chips,
        pheno_list,
        vcf_list,
        f"Pedigree: {n_pedi}",
        f"Phenotype: {n_pheno}",
        f"VCF: {n_vcf}",
    )
    
@callback(
    Output("pedigree-cytoscape", "elements", allow_duplicate=True),
    Output("pedigree-cytoscape", "stylesheet", allow_duplicate=True),
     Output("selected-nodes-store", "data", allow_duplicate=True),  # ✅ 초기화용 Output 추가
    [
        Input("addv-apply-groups-pedi", "data"),
        Input("addv-apply-groups-nopedi", "data"),
    ],
   [
    State("fixed-vcf-item", "data"),
    State("fixed-vcf-item-nopedi", "data"),],
    prevent_initial_call=True,
)
def on_apply_expand_group_pedi_nopedi(group_pedi, group_nopedi, fixed_vcf_item, fixed_vcf_item_nopedi):
    
    """
    - addv-apply-groups-pedi → ±2세대 확장
    - addv-apply-groups-nopedi → 독립 노드 추가 (20개 단위 anchor 클러스터)
      → anchor 중심의 circle 배치 전용
    """
    import traceback, unicodedata, re, math
    from dash.exceptions import PreventUpdate

    if not group_pedi and not group_nopedi:
        raise PreventUpdate

    fixed_nodes = []

    if fixed_vcf_item and isinstance(fixed_vcf_item, dict):
        fixed_nodes.append({
            "id": str(fixed_vcf_item.get("processed_name") or ""),
            "vcf_id": str(fixed_vcf_item.get("status") or ""),
        })

    if fixed_vcf_item_nopedi and isinstance(fixed_vcf_item_nopedi, dict):
        fixed_nodes.append({
            "id": str(fixed_vcf_item_nopedi.get("variety_id") or ""),
            "vcf_id": str(fixed_vcf_item_nopedi.get("status") or ""),
        })

    # ✅ basenode 결정 (둘 중 하나만 사용됨)
    basenode = None
    if fixed_nodes:
        basenode = next((n.get("id") for n in fixed_nodes if n.get("id")), None)
    
    print(f"group_pedi: {group_pedi}")
    print(f"group_nopedi: {group_nopedi}")
    try:
        pedigree_app = get_pedigree_app()
        if not pedigree_app:
            raise PreventUpdate
        merged_elements = []
        fixed_name = None
        if fixed_vcf_item and isinstance(fixed_vcf_item, dict):
            fixed_name = fixed_vcf_item.get("processed_name") or fixed_vcf_item.get("variety_id")

        if fixed_name:
            base_nodes, base_edges = pedigree_app.get_connected_nodes(fixed_name, 2, 2)
            base_elements = pedigree_app.create_cytoscape_elements(base_nodes, base_edges)
            merged_elements.extend(base_elements)
            print(f"🟩 Base pedigree built for {fixed_name} ({len(base_nodes)} nodes)")

        # ----------------------------------------------------
        # 2️⃣ nopedi root: 단일 노드 생성
        # ----------------------------------------------------
        elif fixed_vcf_item_nopedi and isinstance(fixed_vcf_item_nopedi, dict):
            print("🌾 NOPEDI case detected — generating standalone node")

            node_id = str(fixed_vcf_item_nopedi.get("variety_id") or "").strip()
            vcf_status = str(fixed_vcf_item_nopedi.get("status") or "No VCF").strip()

            has_vcf = vcf_status not in ("", "No VCF", None)
            color_class = get_color_class(has_vcf, True)

            merged_elements.append({
                "data": {
                    "id": node_id,
                    "label": node_id,
                    "type": "nopedi",
                    "it_number": node_id,
                    "vcf_id": vcf_status,
                    "vcf_status": vcf_status,
                    "has_vcf": has_vcf,
                    "has_it": True,
                    "color_class": color_class,
                },
                "position": {"x": 0, "y": 0},
                "classes": "nopedi-node",
                "style": {
                    "shape": "rectangle",
                    "width": 90,
                    "height": 30,
                    "font-size": "11px",
                    "text-valign": "center",
                    "text-halign": "center",
                    "border-width": 1.5,
                    "shadow-color": "#7f8c8d",
                    "shadow-blur": 10,
                    "shadow-opacity": 0.6,
                    "shadow-offset-x": 2,
                    "shadow-offset-y": 2,
                    "underlay-color": "#bdc3c7",
                    "underlay-opacity": 0.8,
                    "underlay-padding": 3,
                },
            })

        # ---------------------------------
        # 1️⃣ PEDI 그룹 확장
        # ---------------------------------
        if group_pedi:
            for g in group_pedi:
                gn, ge = pedigree_app.get_connected_nodes(g, 2, 2)
                g_elem = pedigree_app.create_cytoscape_elements(gn, ge)
                merged_elements.extend(g_elem)

        # ---------------------------------
        # 2️⃣ NOPEDI 그룹 (4개 단위 성장형 hourglass + side-anchor)
        # ---------------------------------
        if isinstance(group_nopedi, dict) and "records" in group_nopedi:
            nopedi_records = group_nopedi.get("records", [])
            total = len(nopedi_records)
            print(f"➕ Adding {total} NOPEDI nodes (growth-type hourglass layout)")
            existing_node_ids = {e["data"]["id"] for e in merged_elements if "data" in e and "id" in e["data"]}
            # 안전한 ID 선택
            def _safe_key(rec):
                name = rec.get("name")
                it = rec.get("it_number")
                name = name if isinstance(name, str) else ("" if name is None else str(name))
                it = it if isinstance(it, str) else ("" if it is None else str(it))
                name, it = name.strip(), it.strip()
                if name.lower() in ("none", "null", "nan", ""):
                    name = ""
                if it.lower() in ("none", "null", "nan", ""):
                    it = ""
                return name or it

            # --- 패턴 분류 ---
            batch_size = 4
            patterns = []
            pending_upper = []
            for i in range(0, total, batch_size):
                subset = nopedi_records[i:i + batch_size]
                if not pending_upper:
                    # 첫 4개는 side anchor (독립)
                    pending_upper = subset
                    patterns.append(("side", subset))
                else:
                    # 4개가 이미 있고 또 4개가 오면 → hourglass
                    patterns.append(("hourglass", (pending_upper, subset)))
                    pending_upper = []

            # --- 배치 기본 파라미터 ---
            base_x = 0
            base_y = 1000
            x_spacing = 100
            y_spacing = 150
            side_offset = 250
            anchor_gap = 300
            spread_pattern = [10, 12, 14, 12, 10]
            anchor_index = 1
            anchor_y = base_y
            existing_anchor = None

            # --- 실제 배치 루프 ---
            for ptype, content in patterns:
                if ptype == "side":
                    subset = content
                    anchor_id = f"nopedi_anchor_{anchor_index}"
                    anchor_index += 1

                    # side anchor (오른쪽으로 살짝 이동)
                    merged_elements.append({
                        "data": {"id": anchor_id, "label": f"NOPEDI SIDE-{anchor_index}"},
                        "classes": "nopedi-anchor",
                        "position": {"x": base_x + side_offset, "y": anchor_y},
                        "style": {
                            "background-color": "#ffffff",
                            "border-color": "#555",
                            "border-style": "dashed",
                            "border-width": 2,
                            "shape": "ellipse",
                            "width": 80,
                            "height": 80,
                            "opacity": 0,
                            #"label": f"SIDE-{anchor_index}",
                            "font-size": "10px",
                        },
                    })

                    # 노드 배치 (오른쪽 확장)
                    for i, rec in enumerate(subset):
                        node_key = _safe_key(rec)
                        it_number = rec.get("it_number", "")
                        has_vcf = rec.get("has_vcf", False)
                        has_pheno = rec.get("has_pheno", False)
                        vcf_id = rec.get("vcf_id", "No VCF")
                        color_class = get_color_class(has_vcf, has_pheno)
                        

                        if not node_key or node_key in existing_node_ids:
                            continue
                        x = base_x + side_offset + (i - len(subset)/2) * x_spacing
                        y = anchor_y
                        merged_elements.append({
                            "data": {"id": node_key, "label": node_key, "type": "nopedi", 
                            "it_number": it_number,'vcf_id': vcf_id,'vcf_status': vcf_id ,
                            "has_vcf": has_vcf, "has_it": has_pheno, "color_class": color_class},
                            "position": {"x": x, "y": y},
                            "classes": "nopedi-node",
                            "style": {
                                #"background-color": "#f39c12",
                                "shape": "rectangle",
                                "font-size": "10px",
                                "text-valign": "center",
                                "text-halign": "center",
                            },
                        })
                        merged_elements.append({
                            "data": {"source": anchor_id, "target": node_key},
                            "classes": "nopedi-link",
                            "style": {"line-style": "dotted", "width": 1, "line-color": "#ccc", "opacity": 0},
                        })
                        existing_node_ids.add(node_key)

                    # 다음 anchor 수직 이동
                    anchor_y += anchor_gap

                elif ptype == "hourglass":
                    upper, lower = content
                    anchor_id = f"nopedi_anchor_{anchor_index}"
                    anchor_index += 1

                    # hourglass anchor
                    merged_elements.append({
                        "data": {"id": anchor_id, "label": f"NOPEDI HG-{anchor_index}"},
                        "classes": "nopedi-anchor",
                        "position": {"x": base_x, "y": anchor_y},
                        "style": {
                            "background-color": "#ffffff",
                            "border-color": "#555",
                            "border-style": "dashed",
                            "border-width": 2,
                            "shape": "ellipse",
                            "width": 90,
                            "height": 90,
                            "opacity": 0,
                            "label": f"HG-{anchor_index}",
                            "font-size": "10px",
                        },
                    })

                    # 위쪽 노드
                    for i, rec in enumerate(upper):
                        node_key = _safe_key(rec)
                        it_number = rec.get("it_number", "")
                        has_vcf = rec.get("has_vcf", False)
                        has_pheno = rec.get("has_pheno", False)
                        vcf_id = rec.get("vcf_id", "No VCF")
                        color_class = get_color_class(has_vcf, has_pheno)
                        if not node_key or node_key in existing_node_ids:
                            continue
                        x = base_x + (i - len(upper)/2) * x_spacing
                        y = anchor_y - y_spacing
                        merged_elements.append({
                            "data": {"id": node_key, "label": node_key, "type": "nopedi", 
                            "it_number": it_number,'vcf_id': vcf_id,'vcf_status': vcf_id ,   
                            "has_vcf": has_vcf, "has_it": has_pheno, "color_class": color_class},
                            "position": {"x": x, "y": y},
                            "classes": "nopedi-node",
                            "style": {
                                #"background-color": "#f39c12",
                                "shape": "rectangle",
                                "font-size": "10px",
                                "text-valign": "center",
                                "text-halign": "center",
                            },
                        })
                        merged_elements.append({
                            "data": {"source": node_key, "target": anchor_id},
                            "classes": "nopedi-link",
                            "style": {"line-style": "dotted", "width": 1, "line-color": "#ccc", "opacity": 0},
                        })
                        existing_node_ids.add(node_key)

                    # 아래쪽 노드
                    for i, rec in enumerate(lower):
                        node_key = _safe_key(rec)
                        it_number = rec.get("it_number", "")
                        has_vcf = rec.get("has_vcf", False)
                        has_pheno = rec.get("has_pheno", False)
                        vcf_id = rec.get("vcf_id", "No VCF")
                             
                        color_class = get_color_class(has_vcf, has_pheno)
                        
                        if not node_key or node_key in existing_node_ids:
                            continue
                        x = base_x + (i - len(lower)/2) * x_spacing
                        y = anchor_y + y_spacing
                        merged_elements.append({
                            "data": {"id": node_key, "label": node_key, "type": "nopedi", 
                            "it_number": it_number,'vcf_id': vcf_id,'vcf_status': vcf_id ,    
                            "has_vcf": has_vcf, "has_it": has_pheno, "color_class": color_class},
                            "position": {"x": x, "y": y},
                            "classes": "nopedi-node",
                            "style": {
                                #"background-color": "#f39c12",
                                "shape": "rectangle",
                                "font-size": "10px",
                                "text-valign": "center",
                                "text-halign": "center",
                            },
                        })
                        merged_elements.append({
                            "data": {"source": anchor_id, "target": node_key},
                            "classes": "nopedi-link",
                            "style": {"line-style": "dotted", "width": 1, "line-color": "#ccc", "opacity": 0},
                        })
                        existing_node_ids.add(node_key)

                    # 이전 anchor 아래층 → 현재 anchor 연결 (hourglass 연결)
                    if existing_anchor:
                        merged_elements.append({
                            "data": {"source": existing_anchor, "target": anchor_id},
                            "classes": "anchor-bridge",
                            "style": {"line-style": "dashed", "width": 2, "line-color": "#999", "opacity": 0},
                        })

                    existing_anchor = anchor_id
                    anchor_y += anchor_gap

        print(f"✅ merged {len(merged_elements)} elements total (pedi + nopedi)")

        style=get_default_stylesheet2(basenode)
        cleared_selection = []

        return merged_elements,style,cleared_selection

    except Exception as e:
        print(f"❌ expand error: {e}\n{traceback.format_exc()}")
        raise PreventUpdate


#==detect
def detect_expand_status(elements, path_store, child_store):
    """
    elements를 순회하면서 어떤 노드가 expand / expand-child에 해당하는지 탐지.
    반환: {'expand': set([...]), 'expand_child': set([...])}
    """
    if not elements:
        return {"expand": set(), "expand_child": set()}

    expand_names = {p["name"] for p in (path_store or [])}
    expand_child_names = set()
    for ch_list in (child_store or {}).values():
        expand_child_names.update(ch_list or [])

    element_ids = {
        elem["data"]["id"]
        for elem in elements
        if "data" in elem and "id" in elem["data"]
    }

    expand = expand_names & element_ids
    expand_child = expand_child_names & element_ids

    return {"expand": expand, "expand_child": expand_child}

'''
@callback(
    Output("pedigree-cytoscape", "stylesheet", allow_duplicate=True),
    [
        Input("pedigree-cytoscape", "elements"),
        Input("fixed-vcf-item", "data"),
        Input("fixed-vcf-item-nopedi", "data"),
        Input("available-nodes-store", "data"),
        Input("reset-view-trigger", "data"),   # ✅ 추가
    ],
    [
        State("pedigree-path-store", "data"),
        State("pedigree-path-child-store", "data"),
        State("addv-apply-groups-pedi", "data"),
    ],
    prevent_initial_call=True,
)
def update_styles(elements, fixed_vcf, fixed_vcf_nopedi, available, reset_trigger,
                  path_store, child_store, group_pedi):
    """
    Reset 버튼 클릭 시에도 스타일 재적용되도록 reset_trigger를 감시.
    """
    import unicodedata
    from dash.exceptions import PreventUpdate

    if not elements:
        raise PreventUpdate

    print(f"🎯 update_styles triggered by reset-view-trigger={reset_trigger}")

    base_name = None
    fixed_vcf_data = _extract_fixed_vcf(fixed_vcf, fixed_vcf_nopedi)
    if fixed_vcf_data and isinstance(fixed_vcf, dict):
        base_name = fixed_vcf.get("processed_name") or fixed_vcf.get("variety_id")

    expanded_ids = [p["name"] for p in (path_store or []) if p.get("status") == "active"]
    expanded_child_ids = list({c for _, children in (child_store or {}).items() for c in children})
    available_ids = [a["id"] for a in (available or [])]

    excluded_nodes = set()

    try:
        pedigree_app = get_pedigree_app()
        if pedigree_app is None:
            print("⚠️ pedigree_app is None → skip exclusion")
        else:
            # Base lineage
            if base_name and not is_nopedi_case(fixed_vcf, fixed_vcf_nopedi):
                nodes, edges = pedigree_app.get_connected_nodes(base_name, 2, 2)
                elems = pedigree_app.create_cytoscape_elements(nodes, edges)
                ids = {e["data"]["id"] for e in elems if "source" not in e["data"]}
                excluded_nodes.update(ids)
                print(f"🟩 Base lineage excluded: {len(ids)} nodes")

            # Group lineage
            if group_pedi:
                for g in group_pedi:
                    gkey = unicodedata.normalize("NFC", g.strip().lower().replace("-", "").replace(" ", ""))
                    for n in pedigree_app.graph.nodes:
                        norm = unicodedata.normalize("NFC", n.strip().lower().replace("-", "").replace(" ", ""))
                        if norm == gkey:
                            gn, ge = pedigree_app.get_connected_nodes(n, 2, 2)
                            gelem = pedigree_app.create_cytoscape_elements(gn, ge)
                            gids = {e["data"]["id"] for e in gelem if "source" not in e["data"]}
                            excluded_nodes.update(gids)
                            print(f"🟦 Group '{n}' lineage excluded: {len(gids)} nodes")
                            break

    except Exception as e:
        print(f"❌ lineage node scan error: {e}")

    filtered_child_ids = [nid for nid in expanded_child_ids if nid not in excluded_nodes]

    add_group_labels = {}
    if group_pedi:
        for i, g in enumerate(group_pedi, start=1):
            add_group_labels[g] = f"add item{i}"
    if base_name:
        add_group_labels[base_name] = "default"

    ss = build_stylesheet_with_available2(
        base=[],
        available_ids=available_ids,
        expanded_ids=expanded_ids,
        expanded_child_ids=filtered_child_ids,
        base_name=base_name,
        add_group_labels=add_group_labels,
    )

    print(f"🎨 스타일 재적용 완료 (Reset triggered): base={base_name}, excluded={len(excluded_nodes)}")
    return ss
'''


@callback(
    [
        Output('pedigree-main-container', 'style', allow_duplicate=True),
        
        Output('btn-wide', 'children'),
        Output('pedigree-cytoscape', 'style'),  # ✅ Cytoscape 높이 변경
    ],
    Input('btn-wide', 'n_clicks'),
    prevent_initial_call=True
)
def toggle_wide_view(n_clicks):
    """좌측 → 전체화면 확대/복귀 토글"""
    if n_clicks % 2 == 1:
        # 🔳 확대 모드
        wrapper_style = {
            'position': 'fixed',
            'top': '5%',
            'left': '2%',
            'width': '96%',
            'height': '90%',
            'backgroundColor': '#fff',
            'zIndex': 1000,
            'border': '8px solid rgba(52, 152, 219, 0.6)',  # ✅ 파란 투명 프레임
            'boxShadow': '0 4px 16px rgba(0,0,0,0.3)',
            'borderRadius': '8px',
            'padding': '10px',
            'transition': 'all 0.3s ease-in-out'
        }

        cyto_style = {
            'height': '100%',  # ✅ wrapper 높이에 맞춤
            'width': '100%'
        }

        icon = html.I(className="fas fa-compress")  # 축소 아이콘

    else:
        # 🔲 기본 모드
        wrapper_style = {
                    'padding': '15px',
                    'backgroundColor': '#ffffff',
                    'borderRadius': '8px',
                    'border': '1px solid #dee2e6',
                    'boxShadow': '0 2px 4px rgba(0,0,0,0.05)',
                    'position': 'relative',  # sticky에서 relative로 변경하여 tooltip positioning 지원
                   # 'top': '20px',  # 상단에서 20px 떨어진 위치에 고정
                    'backgroundColor': 'white', # 배경 덮어쓰기
                    'zIndex': 10      
        }

        cyto_style = {
            'height': '600px',  # ✅ 원래 크기로 복귀
            'width': '100%'
        }

        icon = html.I(className="fas fa-expand")

    return wrapper_style, icon, cyto_style

@callback(
    [
        Output('pedigree-cytoscape', 'stylesheet', allow_duplicate=True),
        Output('gwas-click-store', 'data', allow_duplicate=True),
        Output('gwas-sample-scatter', 'figure', allow_duplicate=True),
        Output('gwas-selected-badge', 'children', allow_duplicate=True),
        Output('gwas-sample-table', 'data', allow_duplicate=True),
        Output('gwas-sample-table', 'style_data_conditional', allow_duplicate=True),
        Output('gwas-genotype-legend', 'children', allow_duplicate=True)

    ],
    [
        Input('gwas-sample-scatter', 'clickData'),
        Input('gwas-sample-table', 'active_cell'),
        Input('gwas-selected-samples-store', 'data'),
    ],
    [
        State('gwas-sample-table', 'derived_virtual_data'),
        State('gwas-click-store', 'data'),
        State('gwas-combo-store', 'data'),
        State('gwas-selected-samples-store', 'data'),
        State('gwas-sample-scatter', 'figure'),
        State('gwas-sample-scatter-group-data', 'data'),
        State('selected-nodes-store', 'data'),
        State('fixed-vcf-item', 'data'),
        State('fixed-vcf-item-nopedi', 'data'),
        State('gwas-sample-table', 'page_current'),
        State('gwas-sample-table', 'page_size'),
        State('highlighted-components-store', 'data'),
        State('component-info-store', 'data')
    ],
    prevent_initial_call=True
)
def handle_gwas_click(scatter_click, active_cell, sample_trigger,
                      table_data, click_store, combo_store,
                      selected_samples, fig_state, scatter_meta, 
                      selected_nodes, fixed_vcf_item, fixed_vcf_item_nopedi,
                      page_current, page_size, selected_comps, comp_info):
    from dash import ctx, no_update
    from dash.exceptions import PreventUpdate
    import pandas as pd, re, numpy as np, copy
    import plotly.graph_objects as go

    trigger = ctx.triggered_id
    click_store = click_store or []
    combo_store = combo_store or {}
    variant_only = combo_store.get("variant_only", [])
    selected_samples = selected_samples or []
    selected_nodes = selected_nodes or []

    # ✅ base node 병합
    fixed_nodes = []
    if fixed_vcf_item and isinstance(fixed_vcf_item, dict):
        fixed_nodes.append({
            "id": str(fixed_vcf_item.get("processed_name") or ""),
            "vcf_id": str(fixed_vcf_item.get("status") or ""),
        })
    if fixed_vcf_item_nopedi and isinstance(fixed_vcf_item_nopedi, dict):
        fixed_nodes.append({
            "id": str(fixed_vcf_item_nopedi.get("variety_id") or ""),
            "vcf_id": str(fixed_vcf_item_nopedi.get("status") or ""),
        })
    all_nodes = selected_nodes + fixed_nodes

    # -------------------------------
    # 1️⃣ 샘플 변경 시 리셋
    # -------------------------------
    if trigger == 'gwas-selected-samples-store':
        highlight_style = []
        if selected_comps:
            color_palette = [
                "#8e44ad",  "#6c5ce7", "#2d3436", "#1a1a1a",
                "#c0399f", "#be90d4", "#4b0082", "#5c3c92", "#3d3d3d",
                "#1C6BA0", "#2980b9", "#e74c3c", "#c0392b", "#e67e22",
                "#d35400", "#27ae60", "#2ecc71", "#f1c40f", "#f39c12"
            ]
            for i, comp_id in enumerate(selected_comps):
                comp = next((c for c in comp_info if c.get('id') == comp_id), None)
                if not comp:
                    continue
                ccolor = color_palette[i % len(color_palette)]
                for e in comp.get("edges", []):
                    src, tgt = e.get("source"), e.get("target")
                    if src and tgt:
                        highlight_style.append({
                            "selector": f'edge[source="{src}"][target="{tgt}"]',
                            "style": {
                                "line-color": ccolor,
                                "target-arrow-color": ccolor,
                                "width": 4,
                            },
                        })

        highlight_select_styles = []
        for nd in selected_nodes or []:
            nid = nd.get("id") or nd.get("name")
            if not nid:
                continue
            highlight_select_styles.append({
                "selector": f'node[id = "{nid}"]',
                "style": {
                    "border-color": "#ff33cc",   # 짙은 초록 외곽선
                    "border-width": 6,           # 두꺼운 강조선
                    "transition-property": "border-color, border-width",
                    "transition-duration": "0.3s",
                },
            })
        base_style=get_default_stylesheet2()
        new_styles = base_style + highlight_style +highlight_select_styles

        if not fig_state:
            raise PreventUpdate
        return (
            new_styles,
            [],
            fig_state,
            "Selected Variants (0)",
            table_data,
            []
            ,
            []
        )

    # -------------------------------
    # 2️⃣ variant_id 추출
    # -------------------------------
    variant_id = None
    if trigger == 'gwas-sample-scatter' and scatter_click:
        variant_id = str(scatter_click["points"][0]["customdata"][4])
    elif trigger == 'gwas-sample-table' and active_cell:
        row_idx = active_cell.get("row")
        page_current = int(page_current or 0)
        page_size = int(page_size or 10)
        abs_idx = row_idx + page_current * page_size
        if table_data and abs_idx is not None and 0 <= abs_idx < len(table_data):
            variant_id = str(
                table_data[abs_idx].get("Variant ID")
                or table_data[abs_idx].get("Variant_ID")
            )
    else:
        raise PreventUpdate
    if not variant_id:
        raise PreventUpdate

    # -------------------------------
    # 3️⃣ variant 정보 추출
    # -------------------------------
    df_combo = pd.DataFrame(variant_only)
    if df_combo.empty or "Variant ID" not in df_combo.columns:
        raise PreventUpdate
    row = df_combo.loc[df_combo["Variant ID"].astype(str) == variant_id]
    if row.empty:
        raise PreventUpdate
    row = row.iloc[0].to_dict()
    minor = str(row.get("Minor Allele") or row.get("Minor_Allele") or "")
    trait = str(row.get("Trait") or "")
    chr_ = str(row.get("Chromosome", ""))
    pos = str(row.get("Position", ""))
    Subtrait=str(row.get("Subtrait", ""))

    used_samples = (scatter_meta or {}).get("used_samples", [])
    if used_samples:
        selected_samples = [s for s in selected_samples if s in used_samples]
    
    # ⚙️ used_samples 기반으로 필터링
    valid_nodes = []
    for node in all_nodes:
        vcf_id = str(node.get("vcf_id", "")).strip()
        if not vcf_id or vcf_id not in used_samples:
            continue
        valid_nodes.append(node)


    # -------------------------------
    # 4️⃣ GT 매핑
    # -------------------------------
    sample_gt_map = {s: row.get(f"{s}_GT", "") for s in selected_samples}
    node_map = {}
    for n in valid_nodes:
        key = valid_vcf_value(n.get("vcf_id")) or valid_vcf_value(n.get("vcf_status"))
        if key and n.get("id"):
            node_map[str(key)] = str(n["id"])
    node_gt_map = {node_map.get(s): gt for s, gt in sample_gt_map.items() if s in node_map}

    # -------------------------------
    # 5️⃣ click_store 토글
    # -------------------------------
    new_click_store = click_store.copy()
    existing_ids = {v["Variant_ID"] for v in click_store}
    if variant_id in existing_ids:
        new_click_store = [v for v in click_store if v["Variant_ID"] != variant_id]
    else:
        new_variant = {
            "Variant_ID": variant_id,
            "Chromosome": chr_,
            "Position": pos,
            "Trait": trait,
            "Subtrait":Subtrait,
            "Minor_Allele": minor,
            "Samples": selected_samples,
        }
        new_variant.update({f"{s}_GT": gt for s, gt in sample_gt_map.items()})
        new_click_store.append(new_variant)
    prev_set = {v["Variant_ID"] for v in click_store}
    new_set = {v["Variant_ID"] for v in new_click_store}
    if prev_set == new_set:
        raise PreventUpdate
    click_store = new_click_store
    clicked_ids = list(new_set)

    # -------------------------------
    # 6️⃣ Pedigree 색상 업데이트
    # -------------------------------
    def contains_minor(gt, minor):
        if not gt or not minor:
            return False
        return re.search(rf'\b{re.escape(minor)}\b', gt, re.IGNORECASE) is not None

    node_colors = {}
    if not click_store:
        styles = get_default_stylesheet2()

    else:

        for node in valid_nodes:
            nid = str(node["id"])
            vcf_id = str(node.get("vcf_id", ""))
            matched_sample = next(
                (s for s in selected_samples if s == vcf_id or s in vcf_id or vcf_id in s),
                None
            )
            if not matched_sample:
                node_colors[nid] = "#bcbcbc"
                continue

            node_pass_flags = []
            for var in click_store:
                minor_allele = var.get("Minor_Allele", "")
                gt_value = var.get(f"{matched_sample}_GT", "")
                node_pass_flags.append(contains_minor(gt_value, minor_allele))

            if all(node_pass_flags):
                node_colors[nid] = "#9B111E"  # ✅ 모든 variant minor allele 보유
            elif any(node_pass_flags):
                node_colors[nid] = "#f1c40f"  # ⚠️ 일부만 포함
            else:
                node_colors[nid] = "#323232"  # ❌ 없음
        styles = get_default_stylesheet2() + [
        {"selector": f'[id="{nid}"]',
         "style": {"background-color": color, "border-color": color, "border-width": "3px"}}
        for nid, color in node_colors.items()
    ]

    # -------------------------------
    # 7️⃣ Scatter 업데이트
    # -------------------------------
    if trigger in ('gwas-sample-scatter', 'gwas-sample-table') and fig_state:
        fig = go.Figure(copy.deepcopy(fig_state))
        for trace in fig.data:
            if not getattr(trace, "customdata", None):
                continue
            vids = [str(cd[4]) for cd in trace.customdata]
            orig_colors = trace.marker.color
            if not isinstance(orig_colors, (list, tuple, np.ndarray)):
                orig_colors = [orig_colors] * len(vids)
            trace.marker.color = [
                "#e74c3c" if vid in clicked_ids else orig_colors[i]
                for i, vid in enumerate(vids)
            ]
            trace.marker.opacity = [
                1.0 if vid in clicked_ids else (0.2 if clicked_ids else 1.0)
                for vid in vids
            ]
        scatter_fig = fig
    else:
        scatter_fig = no_update

    # -------------------------------
    # 8️⃣ Table 정렬 + 스타일
    # -------------------------------
    def sort_table_by_click_store(table_data, click_store):
        clicked_ids = [v["Variant_ID"] for v in click_store]
        clicked_set = set(clicked_ids)
        return sorted(
            table_data,
            key=lambda r: (clicked_ids.index(str(r.get("Variant ID")))
                        if str(r.get("Variant ID")) in clicked_set else 9999)
        )
    if table_data and click_store:
        sorted_data = sort_table_by_click_store(table_data, click_store)
    else:
        sorted_data = table_data

    style_conditional = [
        {"if": {"filter_query": f'{{Variant ID}} = "{vid}"'},
         "backgroundColor": "rgba(231,76,60,0.1)",
         "border": "1px solid #e74c3c",
         "fontWeight": "600"}
        for vid in clicked_ids
    ]
    badge_label = f"Selected Variants ({len(click_store)})"



    # --- Component edge 하이라이트 ---
    highlight_style = []
    if selected_comps:
        color_palette = [
            "#8e44ad",  "#6c5ce7", "#2d3436", "#1a1a1a",
            "#c0399f", "#be90d4", "#4b0082", "#5c3c92", "#3d3d3d",
            "#1C6BA0", "#2980b9", "#e74c3c", "#c0392b", "#e67e22",
            "#d35400", "#27ae60", "#2ecc71", "#f1c40f", "#f39c12"
        ]
        for i, comp_id in enumerate(selected_comps):
            comp = next((c for c in comp_info if c.get('id') == comp_id), None)
            if not comp:
                continue
            ccolor = color_palette[i % len(color_palette)]
            for e in comp.get("edges", []):
                src, tgt = e.get("source"), e.get("target")
                if src and tgt:
                    highlight_style.append({
                        "selector": f'edge[source="{src}"][target="{tgt}"]',
                        "style": {
                            "line-color": ccolor,
                            "target-arrow-color": ccolor,
                            "width": 4,
                        },
                    })

    highlight_select_styles = []
    for nd in selected_nodes or []:
        nid = nd.get("id") or nd.get("name")
        if not nid:
            continue
        highlight_select_styles.append({
            "selector": f'node[id = "{nid}"]',
            "style": {
                "border-color": "#ff33cc",   # 짙은 초록 외곽선
                "border-width": 6,           # 두꺼운 강조선
                "transition-property": "border-color, border-width",
                "transition-duration": "0.3s",
            },
        })

    new_styles = styles + highlight_style +highlight_select_styles

    # -------------------------------
    # 9️⃣ Genotype Legend 구성
    # -------------------------------
    legend_title = "Genotype Legend"
    if len(click_store) == 1:
        single = click_store[0]
        minor = single.get("Minor_Allele", "")
        legend_title += f" — Single Variant"
    else:
        legend_title += f" — Multi Variants"

    legend_children = [

        html.Div([
            html.Span("●", style={"color": "#9B111E", "marginRight": "4px"}), "All minor alleles present",
        ], style={"fontSize": "12px", "marginBottom": "2px"}),

        html.Div([
            html.Span("●", style={"color": "#f1c40f", "marginRight": "4px"}), "Partially matched",
        ], style={"fontSize": "12px", "marginBottom": "2px"}),

        html.Div([
            html.Span("●", style={"color": "#323232", "marginRight": "4px"}), "No minor allele",
        ], style={"fontSize": "12px"}),
    ]

    return new_styles, click_store, scatter_fig, badge_label, sorted_data, style_conditional,legend_children

@callback(
    [
        Output('pedigree-cytoscape', 'stylesheet', allow_duplicate=True),
        Output('gwas-click-store', 'data', allow_duplicate=True),
        Output('gwas-sample-scatter-group', 'figure', allow_duplicate=True),
        Output('gwas-selected-badge', 'children', allow_duplicate=True),
        Output('gwas-sample-table-group', 'data', allow_duplicate=True),
        Output('gwas-sample-table-group', 'style_data_conditional', allow_duplicate=True),
        Output('gwas-genotype-legend', 'children', allow_duplicate=True),
    ],
    [
        Input('gwas-sample-scatter-group', 'clickData'),
        Input('gwas-sample-table-group', 'active_cell'),
        Input('gwas-selected-samples-store', 'data'),
    ],
    [
        State('gwas-sample-table-group', 'derived_virtual_data'),
        State('gwas-click-store', 'data'),
        State('gwas-combo-store', 'data'),
        State('gwas-sample-scatter-group', 'figure'),
        State('gwas-sample-scatter-group-data', 'data'),
        State('selected-nodes-store', 'data'),
        State('fixed-vcf-item', 'data'),
        State('fixed-vcf-item-nopedi', 'data'),
        State('gwas-sample-table-group', 'page_current'),
        State('gwas-sample-table-group', 'page_size'),
        State('highlighted-components-store', 'data'),
        State('component-info-store', 'data')
    ],
    prevent_initial_call=True
)
def handle_gwas_click_bygroup(
    scatter_click, active_cell, selected_samples,
    table_data, click_store, combo_store,
    fig_state, scatter_meta, selected_nodes,
    fixed_vcf_item, fixed_vcf_item_nopedi,
    page_current, page_size, selected_comps, comp_info
):
    from dash import ctx, no_update
    from dash.exceptions import PreventUpdate
    import pandas as pd, re, numpy as np, copy
    import plotly.graph_objects as go

    trigger = ctx.triggered_id
    click_store = click_store or []
    combo_store = combo_store or {}
    variant_only = combo_store.get("variant_only", [])
    selected_samples = selected_samples or []
    selected_nodes = selected_nodes or []

    # -------------------------------
    # (0) base node 병합
    # -------------------------------
    fixed_nodes = []
    if fixed_vcf_item and isinstance(fixed_vcf_item, dict):
        fixed_nodes.append({
            "id": str(fixed_vcf_item.get("processed_name") or ""),
            "vcf_id": str(fixed_vcf_item.get("status") or ""),
        })
    if fixed_vcf_item_nopedi and isinstance(fixed_vcf_item_nopedi, dict):
        fixed_nodes.append({
            "id": str(fixed_vcf_item_nopedi.get("variety_id") or ""),
            "vcf_id": str(fixed_vcf_item_nopedi.get("status") or ""),
        })
    all_nodes = selected_nodes + fixed_nodes

    # -------------------------------
    # (1) 샘플 변경 → 초기화
    # -------------------------------
    if trigger == 'gwas-selected-samples-store':
        highlight_style = []
        if selected_comps:
            color_palette = [
                "#8e44ad",  "#6c5ce7", "#2d3436", "#1a1a1a",
                "#c0399f", "#be90d4", "#4b0082", "#5c3c92", "#3d3d3d",
                "#1C6BA0", "#2980b9", "#e74c3c", "#c0392b", "#e67e22",
                "#d35400", "#27ae60", "#2ecc71", "#f1c40f", "#f39c12"
            ]
            for i, comp_id in enumerate(selected_comps):
                comp = next((c for c in comp_info if c.get('id') == comp_id), None)
                if not comp:
                    continue
                ccolor = color_palette[i % len(color_palette)]
                for e in comp.get("edges", []):
                    src, tgt = e.get("source"), e.get("target")
                    if src and tgt:
                        highlight_style.append({
                            "selector": f'edge[source="{src}"][target="{tgt}"]',
                            "style": {
                                "line-color": ccolor,
                                "target-arrow-color": ccolor,
                                "width": 4,
                            },
                        })

        highlight_select_styles = []
        for nd in selected_nodes or []:
            nid = nd.get("id") or nd.get("name")
            if not nid:
                continue
            highlight_select_styles.append({
                "selector": f'node[id = "{nid}"]',
                "style": {
                    "border-color": "#ff33cc",   # 짙은 초록 외곽선
                    "border-width": 6,           # 두꺼운 강조선
                    "transition-property": "border-color, border-width",
                    "transition-duration": "0.3s",
                },
            })
        base_style=get_default_stylesheet2()
        new_styles = base_style + highlight_style +highlight_select_styles
        if not fig_state:
            raise PreventUpdate
        return (
            new_styles,
            [],
            fig_state,
            "Selected Variants (0)",
            table_data,
            [],
            []
        )

    # -------------------------------
    # (2) variant_id 추출
    # -------------------------------
    variant_id = None
    if trigger == 'gwas-sample-scatter-group' and scatter_click:
        variant_id = str(scatter_click["points"][0]["customdata"][4])
    elif trigger == 'gwas-sample-table-group' and active_cell:
        row_idx = active_cell.get("row")
        page_current = int(page_current or 0)
        page_size = int(page_size or 10)
        abs_idx = row_idx + page_current * page_size
        if table_data and abs_idx is not None and 0 <= abs_idx < len(table_data):
            variant_id = str(
                table_data[abs_idx].get("Variant ID")
                or table_data[abs_idx].get("Variant_ID")
            )
    if not variant_id:
        raise PreventUpdate

    # -------------------------------
    # (3) variant 정보 lookup
    # -------------------------------
    df_combo = pd.DataFrame(variant_only)
    if df_combo.empty or "Variant ID" not in df_combo.columns:
        raise PreventUpdate
    row = df_combo.loc[df_combo["Variant ID"].astype(str) == variant_id]
    if row.empty:
        raise PreventUpdate
    row = row.iloc[0].to_dict()

    minor = str(row.get("Minor Allele") or "")
    trait = str(row.get("Trait") or "")
    chr_ = str(row.get("Chromosome", ""))
    pos = str(row.get("Position", ""))
    Subtrait=str(row.get("Subtrait", ""))

    # -------------------------------
    # (4) sample GT 매핑
    # -------------------------------
    used_samples = (scatter_meta or {}).get("used_samples", [])
    if used_samples:
        selected_samples = [s for s in selected_samples if s in used_samples]

    valid_nodes = []
    for node in all_nodes:
        vcf_id = str(node.get("vcf_id", "")).strip()
        if not vcf_id or vcf_id not in used_samples:
            continue
        valid_nodes.append(node)

    sample_gt_map = {s: row.get(f"{s}_GT", "") for s in selected_samples}
    node_map = {}
    for n in valid_nodes:
        key = valid_vcf_value(n.get("vcf_id")) or valid_vcf_value(n.get("vcf_status"))
        if key and n.get("id"):
            node_map[str(key)] = str(n["id"])
    node_gt_map = {node_map.get(s): gt for s, gt in sample_gt_map.items() if s in node_map}

    # -------------------------------
    # (5) click_store toggle
    # -------------------------------
    existing_ids = {v["Variant_ID"] for v in click_store}
    if variant_id in existing_ids:
        click_store = [v for v in click_store if v["Variant_ID"] != variant_id]
    else:
        new_entry = {
            "Variant_ID": variant_id,
            "Chromosome": chr_,
            "Position": pos,
            "Subtrait":Subtrait,
            "Trait": trait,
            "Minor_Allele": minor,
            "Samples": selected_samples,
        }
        new_entry.update({f"{s}_GT": gt for s, gt in sample_gt_map.items()})
        click_store.append(new_entry)

    clicked_ids = [v["Variant_ID"] for v in click_store]

    # -------------------------------
    # (6) Pedigree 색상 계산 (GT 기반)
    # -------------------------------
    def contains_minor(gt, minor):
        if not gt or not minor:
            return False
        return re.search(rf'\b{re.escape(minor)}\b', gt, re.IGNORECASE) is not None

    node_colors = {}
    if not click_store:
        styles = get_default_stylesheet2()

    else:

        for node in valid_nodes:
            nid = str(node["id"])
            vcf_id = str(node.get("vcf_id", ""))
            matched_sample = next(
                (s for s in selected_samples if s == vcf_id or s in vcf_id or vcf_id in s),
                None
            )
            if not matched_sample:
                node_colors[nid] = "#bcbcbc"
                continue

            node_pass_flags = []
            for var in click_store:
                minor_allele = var.get("Minor_Allele", "")
                gt_value = var.get(f"{matched_sample}_GT", "")
                node_pass_flags.append(contains_minor(gt_value, minor_allele))

            if all(node_pass_flags):
                node_colors[nid] = "#9B111E"  # ✅ 모든 variant minor allele 보유
            elif any(node_pass_flags):
                node_colors[nid] = "#f1c40f"  # ⚠️ 일부만 포함
            else:
                node_colors[nid] = "#323232"  # ❌ 없음
        styles = get_default_stylesheet2() + [
            {"selector": f'[id="{nid}"]',
            "style": {"background-color": color, "border-color": color, "border-width": "3px"}}
            for nid, color in node_colors.items()
        ]

    # -------------------------------
    # (7) scatter 업데이트
    # -------------------------------
    scatter_fig = no_update
    if fig_state and fig_state.get("data"):
        fig = go.Figure(copy.deepcopy(fig_state))
        for trace in fig.data:
            if not getattr(trace, "customdata", None):
                continue
            vids = [str(cd[4]) for cd in trace.customdata]
            orig_colors = trace.marker.color
            if not isinstance(orig_colors, (list, tuple, np.ndarray)):
                orig_colors = [orig_colors] * len(vids)
            trace.marker.color = [
                "#e74c3c" if vid in clicked_ids else orig_colors[i]
                for i, vid in enumerate(vids)
            ]
            trace.marker.opacity = [
                1.0 if vid in clicked_ids else (0.2 if clicked_ids else 1.0)
                for vid in vids
            ]
        scatter_fig = fig

    # -------------------------------
    # (8) Table 정렬 + 스타일
    # -------------------------------
    def sort_table_by_click_store(table_data, click_store):
        clicked_ids = [v["Variant_ID"] for v in click_store]
        return sorted(
            table_data,
            key=lambda r: (clicked_ids.index(str(r.get("Variant ID")))
                           if str(r.get("Variant ID")) in clicked_ids else 9999)
        )

    sorted_data = sort_table_by_click_store(table_data, click_store) if table_data else table_data
    style_conditional = [
        {"if": {"filter_query": f'{{Variant ID}} = "{vid}"'},
         "backgroundColor": "rgba(231,76,60,0.1)",
         "border": "1px solid #e74c3c",
         "fontWeight": "600"}
        for vid in clicked_ids
    ]
    badge_label = f"Selected Variants ({len(click_store)})"
    '''
    styles = get_default_stylesheet2() + [
        {"selector": f'[id="{nid}"]',
         "style": {"background-color": color, "border-color": color, "border-width": "3px"}}
        for nid, color in node_colors.items()
    ]'''

    # --- Component edge 하이라이트 ---
    highlight_style = []
    if selected_comps:
        color_palette = [
            "#8e44ad",  "#6c5ce7", "#2d3436", "#1a1a1a",
            "#c0399f", "#be90d4", "#4b0082", "#5c3c92", "#3d3d3d",
            "#1C6BA0", "#2980b9", "#e74c3c", "#c0392b", "#e67e22",
            "#d35400", "#27ae60", "#2ecc71", "#f1c40f", "#f39c12"
        ]
        for i, comp_id in enumerate(selected_comps):
            comp = next((c for c in comp_info if c.get('id') == comp_id), None)
            if not comp:
                continue
            ccolor = color_palette[i % len(color_palette)]
            for e in comp.get("edges", []):
                src, tgt = e.get("source"), e.get("target")
                if src and tgt:
                    highlight_style.append({
                        "selector": f'edge[source="{src}"][target="{tgt}"]',
                        "style": {
                            "line-color": ccolor,
                            "target-arrow-color": ccolor,
                            "width": 4,
                        },
                    })
    highlight_select_styles = []
    for nd in selected_nodes or []:
        nid = nd.get("id") or nd.get("name")
        if not nid:
            continue
        highlight_select_styles.append({
            "selector": f'node[id = "{nid}"]',
            "style": {
                "border-color": "#ff33cc",   # 짙은 초록 외곽선
                "border-width": 6,           # 두꺼운 강조선
                "transition-property": "border-color, border-width",
                "transition-duration": "0.3s",
            },
        })

    new_styles = styles + highlight_style +highlight_select_styles

    # -------------------------------
    # 9️⃣ Genotype Legend 구성
    # -------------------------------
    legend_title = "Genotype Legend"
    if len(click_store) == 1:
        single = click_store[0]
        minor = single.get("Minor_Allele", "")
        legend_title += f" — Single Variant ({minor})"
    else:
        legend_title += f" — Multi Variants ({len(click_store)})"

    legend_children = [


        html.Div([
            html.Span("●", style={"color": "#9B111E", "marginRight": "4px"}), "All minor alleles present",
        ], style={"fontSize": "12px", "marginBottom": "2px"}),

        html.Div([
            html.Span("●", style={"color": "#f1c40f", "marginRight": "4px"}), "Partially matched",
        ], style={"fontSize": "12px", "marginBottom": "2px"}),

        html.Div([
            html.Span("●", style={"color": "#323232", "marginRight": "4px"}), "No minor allele",
        ], style={"fontSize": "12px"}),
    ]

    return new_styles, click_store, scatter_fig, badge_label, sorted_data, style_conditional ,legend_children

'''
@callback(
    Output('legend-container', 'children'),
    Input('gwas-click-store', 'data'),
    prevent_initial_call=False
)
def update_legend_based_on_gwas(click_store):
    from dash import html

    # ✅ GWAS variant 선택됨 → 2x2 grid 스타일 legend
    if click_store and len(click_store) > 0:
        print("🎨 GWAS Legend activated")
        return html.Div([
            html.H6("🧬 GWAS Variant Legend",
                    style={'margin': '0 0 8px 0', 'fontSize': '14px', 'fontWeight': 'bold'}),

            html.Div([
                _legend_circle("#27ae60", "#27ae60", "Contains minor allele"),
                _legend_circle("#7f8c8d", "#7f8c8d", "Has GT (no minor)"),
                _legend_circle("#bcbcbc", "#bcbcbc", "No GT"),
                _legend_circle("#6FB1FC", "#5A9FEB", "Basic Node (default)")
            ],
                style={
                    "display": "grid",
                    "gridTemplateColumns": "repeat(2, auto)",  # ✅ 2열 배치
                    "gap": "6px 12px",
                    "alignItems": "center",
                    "justifyContent": "start"
                }
            )
        ],
            style={
                'position': 'relative',
                'padding': '8px 10px',
                'backgroundColor': 'rgba(255, 255, 255, 0.95)',
                'border': '1px solid #dee2e6',
                'borderRadius': '6px',
                'boxShadow': '0 2px 4px rgba(0,0,0,0.1)',
                'width': '240px',
                'marginBottom': '3px'
            }
        )

    # ✅ 기본 Pedigree Legend 복귀
    print("🎨 Pedigree Legend restored")
    return html.Div([
        html.H6("🏷️ Pedigree Legend",
                style={'margin': '0 0 8px 0', 'fontSize': '14px', 'fontWeight': 'bold'}),

        html.Div([
            _legend_circle("#6FB1FC", "#5A9FEB", "Basic Node"),
            _legend_circle("#e74c3c", "#c0392b", "Available Node (VCF / Phenotype)"),
            _legend_circle("#ffffff", "#27ae60", "Expanded Node (in Pedigree Path)"),
            _legend_circle("#f39c12", "#e67e22", "Expanded Child Node"),
            _legend_circle("#f39c12", "#27ae60", "Expanded + Child Node"),
            _legend_circle("#e74c3c", "#212121", "Default Node (in Pedigree Path, Available)"),
            _legend_circle("#424949", "#212121", "Default Node (in Pedigree Path, No Data)"),
            _legend_circle("#e74c3c", "#212121", "Add Item Node (in Pedigree Path, Available)"),
            _legend_circle("#424949", "#212121", "Add Item Node (in Pedigree Path, No Data)"),
        ])
    ],
        style={
            'position': 'relative',
            'padding': '8px 10px',
            'backgroundColor': 'rgba(255, 255, 255, 0.95)',
            'border': '1px solid #dee2e6',
            'borderRadius': '6px',
            'boxShadow': '0 2px 4px rgba(0,0,0,0.1)',
            'width': '240px',
            'marginBottom': '3px'
        }
    )
'''
#on_tap_수정

@callback(
    Output('selected-nodes-store', 'data'),
    Input('pedigree-cytoscape', 'tapNodeData'),
    State('selected-nodes-store', 'data'),
    prevent_initial_call=True
)
def on_tap_select(nd, selected_nodes):
    """노드 클릭 시 선택/해제 (dict 저장, multi-select 지원)"""
    from dash.exceptions import PreventUpdate
    if not nd or 'id' not in nd:
        raise PreventUpdate

    selected_nodes = selected_nodes or []
    nid = nd['id']

    # 현재 선택된 노드들의 id 리스트 추출
    current_ids = [x.get('id') for x in selected_nodes]

    if nid in current_ids:
        # 이미 있으면 제거
        selected_nodes = [x for x in selected_nodes if x.get('id') != nid]
        print(f"🔹 Unselected: {nid}")
    else:
        # 없으면 추가
        selected_nodes.append(nd)
        print(f"✅ Selected: {nid}")

    return selected_nodes

@callback(
    Output('pedigree-cytoscape', 'stylesheet', allow_duplicate=True),
    [
        Input('selected-nodes-store', 'data'),
    ],
    [
        State('fixed-vcf-item', 'data'),
        State('fixed-vcf-item-nopedi', 'data'),
        State('highlighted-components-store', 'data'),
        State('component-info-store', 'data'),
    ],
    prevent_initial_call=True
)
def highlight_selected_nodes(selected_nodes, fixed_vcf_item, fixed_vcf_item_nopedi,selected_comps, comp_info):
    """
    기본 스타일(get_default_stylesheet2) + 선택 노드 외곽선 강조 (짙은 초록)
    """
    # ✅ 기준 노드(basenode) 자동 결정
    basenode = None
    if fixed_vcf_item and isinstance(fixed_vcf_item, dict):
        basenode = fixed_vcf_item.get("processed_name") or fixed_vcf_item.get("variety_id")
    elif fixed_vcf_item_nopedi and isinstance(fixed_vcf_item_nopedi, dict):
        basenode = fixed_vcf_item_nopedi.get("variety_id")

    # ✅ 1️⃣ 기본 스타일
    base_styles = get_default_stylesheet2(basenode)

    # ✅ 2️⃣ 선택 노드 외곽선 강조
    highlight_styles = []
    for nd in selected_nodes or []:
        nid = nd.get("id") or nd.get("name")
        if not nid:
            continue
        highlight_styles.append({
            "selector": f'node[id = "{nid}"]',
            "style": {
                "border-color": "#ff33cc",   # 짙은 초록 외곽선
                "border-width": 6,           # 두꺼운 강조선
                "transition-property": "border-color, border-width",
                "transition-duration": "0.3s",
            },
        })
    
    highlight_style2 = []
    if selected_comps:
        color_palette = [
            "#8e44ad",  "#6c5ce7", "#2d3436", "#1a1a1a",
            "#c0399f", "#be90d4", "#4b0082", "#5c3c92", "#3d3d3d",
            "#1C6BA0", "#2980b9", "#e74c3c", "#c0392b", "#e67e22",
            "#d35400", "#27ae60", "#2ecc71", "#f1c40f", "#f39c12"
        ]
        for i, comp_id in enumerate(selected_comps):
            comp = next((c for c in comp_info if c.get('id') == comp_id), None)
            if not comp:
                continue
            ccolor = color_palette[i % len(color_palette)]
            
            for e in comp.get("edges", []):
                src, tgt = e.get("source"), e.get("target")
                if src and tgt:
                    highlight_style2.append({
                        "selector": f'edge[source="{src}"][target="{tgt}"]',
                        "style": {
                            "line-color": ccolor,
                            "target-arrow-color": ccolor,
                            "width": 4,
                        },
                    })


    # ✅ 3️⃣ 병합 후 반환
    return base_styles + highlight_styles + highlight_style2


# ----------------------------------------
# 2️⃣ 선택 노드 확장 (expand)
# ----------------------------------------
@callback(
    [
        Output('pedigree-cytoscape', 'elements', allow_duplicate=True),
        Output('pedigree-path-store', 'data', allow_duplicate=True),
        Output('pedigree-path-child-store', 'data', allow_duplicate=True),
        Output('pedigree-cytoscape', 'stylesheet', allow_duplicate=True),
    ],
    Input('btn-expand', 'n_clicks'),
    [
        State('selected-nodes-store', 'data'),
        State('pedigree-cytoscape', 'elements'),
        State('pedigree-path-store', 'data'),
        State('pedigree-path-child-store', 'data'),
        State('highlighted-components-store', 'data'),
        State('component-info-store', 'data')
        ],
    prevent_initial_call=True
)
def expand_selected_nodes(n_clicks, selected_nodes, base_elements, path_store, child_store,
                            selected_comps, comp_info):
    """선택된 노드를 확장 (+/- 1 generation)"""
    from dash.exceptions import PreventUpdate
    if not n_clicks or not selected_nodes:
        raise PreventUpdate

    pedigree_app = get_pedigree_app()
    merged = base_elements.copy()
    path_store = path_store or []
    child_store = child_store or {}

    # ✅ dict 구조 대응
    for nd in selected_nodes:
        nid = nd.get("id")
        if not nid:
            continue
        try:
            new_nodes, new_edges = pedigree_app.get_connected_nodes(nid, 1, 1)
            new_elements = pedigree_app.create_cytoscape_elements(new_nodes, new_edges)

            existing_ids = {e["data"]["id"] for e in merged if "source" not in e["data"]}
            for ne in new_elements:
                d = ne.get("data", {})
                if "source" in d:
                    if not any(
                        e.get("data", {}).get("source") == d["source"]
                        and e.get("data", {}).get("target") == d["target"]
                        for e in merged
                    ):
                        merged.append(ne)
                else:
                    if d.get("id") not in existing_ids:
                        merged.append(ne)

            child_store[nid] = sorted(set(child_store.get(nid, [])) | set(new_nodes) - {nid})
            print(f"🌱 Expanded: {nid}")

        except Exception as e:
            print(f"❌ Expand error for {nid}: {e}")

    base_stylesheet = get_default_stylesheet2()
    highlight_style = []
    if selected_comps:
        color_palette = [
            "#8e44ad",  "#6c5ce7", "#2d3436", "#1a1a1a",
            "#c0399f", "#be90d4", "#4b0082", "#5c3c92", "#3d3d3d",
            "#1C6BA0", "#2980b9", "#e74c3c", "#c0392b", "#e67e22",
            "#d35400", "#27ae60", "#2ecc71", "#f1c40f", "#f39c12"
        ]
        for i, comp_id in enumerate(selected_comps):
            comp = next((c for c in comp_info if c.get('id') == comp_id), None)
            if not comp:
                continue
            ccolor = color_palette[i % len(color_palette)]
            for e in comp.get("edges", []):
                src, tgt = e.get("source"), e.get("target")
                if src and tgt:
                    highlight_style.append({
                        "selector": f'edge[source="{src}"][target="{tgt}"]',
                        "style": {
                            "line-color": ccolor,
                            "target-arrow-color": ccolor,
                            "width": 4,
                        },
                    })
    
    highlight_select_styles = []
    for nd in selected_nodes or []:
        nid = nd.get("id") or nd.get("name")
        if not nid:
            continue
        highlight_select_styles.append({
            "selector": f'node[id = "{nid}"]',
            "style": {
                "border-color": "#ff33cc",   # 짙은 초록 외곽선
                "border-width": 6,           # 두꺼운 강조선
                "transition-property": "border-color, border-width",
                "transition-duration": "0.3s",
            },
        })

    return merged, path_store, child_store, base_stylesheet + highlight_style + highlight_select_styles



# ----------------------------------------
# 3️⃣ 선택 노드 제거 (remove)
# ----------------------------------------
@callback(
    [
        Output('pedigree-cytoscape', 'elements', allow_duplicate=True),
        Output('selected-nodes-store', 'data', allow_duplicate=True),
    ],
    Input('btn-remove', 'n_clicks'),
    [
        State('selected-nodes-store', 'data'),
        State('pedigree-cytoscape', 'elements'),
        State('fixed-vcf-item', 'data'),
        State('fixed-vcf-item-nopedi', 'data'),
    ],
    prevent_initial_call=True
)
def remove_selected_nodes(n_clicks, selected_nodes, elements,
                          fixed_vcf_item, fixed_vcf_item_nopedi):
    """선택된 노드를 Cytoscape 그래프에서 제거하되, basenode는 절대 삭제하지 않음"""
    from dash.exceptions import PreventUpdate
    if not n_clicks or not selected_nodes:
        raise PreventUpdate

    # ✅ 1️⃣ 기준 노드(basenode) 확인
    fixed_nodes = []

    if fixed_vcf_item and isinstance(fixed_vcf_item, dict):
        fixed_nodes.append({
            "id": str(fixed_vcf_item.get("processed_name") or ""),
            "vcf_id": str(fixed_vcf_item.get("status") or "")
        })

    if fixed_vcf_item_nopedi and isinstance(fixed_vcf_item_nopedi, dict):
        fixed_nodes.append({
            "id": str(fixed_vcf_item_nopedi.get("variety_id") or ""),
            "vcf_id": str(fixed_vcf_item_nopedi.get("status") or "")
        })

    basenode = None
    if fixed_nodes:
        basenode = next((n.get("id") for n in fixed_nodes if n.get("id")), None)

    # ✅ 2️⃣ 제거 대상 ID 수집 (basenode는 제외)
    selected_ids = [
        nd.get("id") for nd in (selected_nodes or [])
        if nd.get("id") and nd.get("id") != basenode
    ]

    if not selected_ids:
        print("⚠️ No removable nodes (only base node selected).")
        raise PreventUpdate

    print(f"🗑 Removing {len(selected_ids)} nodes: {selected_ids}")
    if basenode:
        print(f"🔒 Base node '{basenode}' will be preserved.")

    # ✅ 3️⃣ elements에서 선택 노드 및 연결 엣지 제거
    new_elements = [
        e for e in elements
        if not (
            ("data" in e and e["data"].get("id") in selected_ids)
            or ("data" in e and (
                e["data"].get("source") in selected_ids or
                e["data"].get("target") in selected_ids)
            )
        )
    ]

    # ✅ 4️⃣ 선택 초기화
    cleared_selection = [
        
    ]

    print("🧹 Cleared selection store after removal (base node kept).")
    return new_elements, cleared_selection


@callback(
    [
        Output('component-info-store', 'data'),
        Output('component-badge-container', 'children'),
    ],
    [
        Input('pedigree-cytoscape', 'elements'),
        Input('reset-view-trigger', 'data'),   # ✅ 초기 로딩 및 리셋 이벤트 추가
    ],
    prevent_initial_call=False   # ✅ 초기 렌더링도 실행
)
def analyze_connected_components(elements, reset_trigger):
    """elements에서 connected components 계산 후 badge 생성 (nopedi anchor 완전 필터링, id기반 처리 포함)"""
    import networkx as nx
    from dash import html

    # ----------------------------
    # ✅ 0️⃣ elements 존재 확인
    # ----------------------------
    if not elements or len(elements) == 0:
        print("⚠️ No elements yet — skipping component analysis.")
        return [], []

    G = nx.Graph()
    node_types = {}
    anchor_nodes = set()

    # ----------------------------
    # 1️⃣ 1차 스캔 — anchor id 수집
    # ----------------------------
    for e in elements:
        d = e.get("data", {})
        classes = e.get("classes", "")
        node_id = d.get("id")

        if not node_id:
            continue

        # classes로 anchor 탐지
        if "nopedi-anchor" in classes:
            anchor_nodes.add(node_id)
            continue

        # id 기반 anchor 패턴 탐지
        if str(node_id).lower().startswith("nopedi_anchor_"):
            anchor_nodes.add(node_id)

    print(f"🔍 Anchors detected: {len(anchor_nodes)} → {sorted(anchor_nodes)}")

    # ----------------------------
    # 2️⃣ 2차 스캔 — 실제 그래프 구성
    # ----------------------------
    for e in elements:
        d = e.get("data", {})
        node_id = d.get("id")
        classes = e.get("classes", "")

        # anchor 노드는 추가하지 않음
        if node_id in anchor_nodes:
            continue

        # anchor-bridge edge 제외
        if "anchor-bridge" in classes:
            continue

        # edge 추가 (anchor 포함 edge 제외)
        if "source" in d and "target" in d:
            if d["source"] in anchor_nodes or d["target"] in anchor_nodes:
                continue
            G.add_edge(d["source"], d["target"])
        elif "id" in d:
            node_types[node_id] = d.get("type", "pedigree")
            G.add_node(node_id)

    # ----------------------------
    # 3️⃣ 컴포넌트 계산
    # ----------------------------
    components = list(nx.connected_components(G))
    comp_info, badges = [], []
    nopedi_nodes = {n for n, t in node_types.items() if t == "nopedi"}

    comp_index = 1
    for comp in components:
        comp_pure_nodes = [n for n in comp if n not in anchor_nodes and n not in nopedi_nodes]
        if not comp_pure_nodes:
            continue

        sub_edges = [
            e for e in elements
            if "data" in e and "source" in e["data"]
            and e["data"]["source"] in comp_pure_nodes
            and e["data"]["target"] in comp_pure_nodes
        ]
        comp_data = {
            "id": f"comp-{comp_index}",
            "nodes": list(comp_pure_nodes),
            "edges": [e["data"] for e in sub_edges],
            "type": "pedigree",
        }
        comp_info.append(comp_data)

        badges.append(
            html.Span(
                f"Pedigree {comp_index} : {len(comp_pure_nodes)} samples",
                id={'type': 'component-badge', 'index': f"comp-{comp_index}"},
                n_clicks=0,
                style={
                    "display": "inline-block",
                            "padding": "4px 10px",
                            "margin": "3px 6px 3px 0",
                            "borderRadius": "6px",
                            "fontSize": "12px",
                            "cursor": "pointer",
                            "transition": "all 0.2s ease-in-out",
                            "userSelect": "none",
                            "whiteSpace": "nowrap",
                            "background": "#9ca3af",   # ⚙️ 기본 회색
                            "color": "#333",
                            "fontWeight": "400",
                            "boxShadow": "none",
                            "transform": "scale(1.0)",
                }
            )
        )
        comp_index += 1

    # ----------------------------
    # 4️⃣ NOPEDI 카운트
    # ----------------------------
    nopedi_count = len(nopedi_nodes)
    if nopedi_count > 0:
        comp_info.append({
            "id": "comp-nopedi",
            "nodes": list(nopedi_nodes),
            "edges": [],
            "type": "nopedi",
        })

        badges.append(
            html.Span(
                f"NOPEDI : {nopedi_count} samples",
                id={'type': 'component-badge', 'index': "comp-nopedi"},
                n_clicks=0,
                style={
                    "display": "inline-block",
        "padding": "4px 10px",
        "margin": "3px 6px 3px 0",
        "borderRadius": "6px",
        "fontSize": "12px",            
        "cursor": "default",          # ✅ 클릭커서 제거
        "pointerEvents": "none",      # ✅ 클릭 이벤트 완전 차단
        "transition": "all 0.2s ease-in-out",
        "userSelect": "none",
        "whiteSpace": "nowrap",
        "background": "#f4b860",  
        "color": "#3b2600",
        "fontWeight": "400",
        "boxShadow": "none",
        "transform": "scale(1.0)",
                }
            )
        )

    print(f"🧩 Found {len(comp_info)} components (Pedigree + NOPEDI), anchors skipped {len(anchor_nodes)})")
    return comp_info, badges



# ---------------------------
# Highlight callback
# ---------------------------
@dash.callback(
    [
        Output('highlighted-components-store', 'data'),
        Output({'type': 'component-badge', 'index': ALL}, 'style'),
    ],
    Input({'type': 'component-badge', 'index': ALL}, 'n_clicks'),
    State('highlighted-components-store', 'data'),
    State('component-info-store', 'data'),
    prevent_initial_call=True
)
def toggle_component_highlight(clicks, highlighted, comp_info):
    from dash import ctx
    from dash.exceptions import PreventUpdate

    if not clicks or not comp_info:
        raise PreventUpdate

    highlighted = highlighted or []
    triggered = [i for i, c in enumerate(clicks) if c and c > 0]
    if not triggered:
        raise PreventUpdate

    comp_id = eval(ctx.triggered[0]['prop_id'].split('.')[0])['index']

    # ✅ 토글
    if comp_id in highlighted:
        highlighted.remove(comp_id)
        print(f"🔹 Unhighlighted {comp_id}")
    else:
        highlighted.append(comp_id)
        print(f"✅ Highlighted {comp_id}")

    # ✅ 스타일 구성
    styles = []
    for cinfo in comp_info:
        cid = cinfo.get("id")
        is_active = cid in highlighted

        base_style = {
            "display": "inline-block",
            "padding": "4px 10px",
            "margin": "3px 6px 3px 0",
            "borderRadius": "6px",
            "fontSize": "12px",
            "cursor": "pointer",
            "transition": "all 0.2s ease-in-out",
            "userSelect": "none",
            "whiteSpace": "nowrap",
        }

        if is_active:
            base_style.update({
                "background": "#1C6BA0",
                "color": "white",
                "fontWeight": "600",
                "boxShadow": "0 2px 5px rgba(0,0,0,0.25)",
                "transform": "scale(1.03)",
            })
        else:
            base_style.update({
                "background": "#9ca3af",  # 밝은 회색
                "color": "#333",
                "fontWeight": "400",
                "boxShadow": "none",
                "transform": "scale(1.0)",
            })
        styles.append(base_style)

    return highlighted, styles


@callback(
    Output('pedigree-cytoscape', 'stylesheet', allow_duplicate=True),
    [Input('highlighted-components-store', 'data')],
    [State('component-info-store', 'data')],
    prevent_initial_call=True
)
def highlight_components(selected_comps, comp_info):
    """
    기본 color_class 색상은 유지하고,
    선택된 컴포넌트에는 외곽선/그림자/underlay 효과만 추가.
    """
    from dash.exceptions import PreventUpdate

    # 1️⃣ 기본 스타일 (데이터 유무 색상 포함)
    base_style = get_default_stylesheet2()

    if not selected_comps:
        return base_style

    # 2️⃣ 하이라이트용 color palette
    color_palette = [
            "#8e44ad",  "#6c5ce7", "#2d3436", "#1a1a1a",
            "#c0399f", "#be90d4", "#4b0082", "#5c3c92", "#3d3d3d",
            "#1C6BA0", "#2980b9", "#e74c3c", "#c0392b", "#e67e22",
            "#d35400", "#27ae60", "#2ecc71", "#f1c40f", "#f39c12"
        ]

    highlight_style = []
    #print(selected_comps)
    # 3️⃣ 선택된 컴포넌트마다 하이라이트 레이어 추가
    for i, comp_id in enumerate(selected_comps):
        comp = next((c for c in comp_info if c.get('id') == comp_id), None)
        if not comp:
            continue

        ccolor = color_palette[i % len(color_palette)]
        #print(ccolor)
        #print(i)
        #print(comp_id)

        #print(comp)
        # 🌟 엣지 하이라이트
        for e in comp.get("edges", []):
            src, tgt = e.get("source"), e.get("target")
            if not src or not tgt:
                continue
            highlight_style.append({
                "selector": f'edge[source = "{src}"][target = "{tgt}"]',
                "style": {
                    "line-color": ccolor,
                    "target-arrow-color": ccolor,
                    "width": 4,
                },
            })

    # ✅ 기본 스타일 + 하이라이트 레이어 결합
    return base_style + highlight_style

@callback(
    [
        Output('gwas-selected-samples-store', 'data', allow_duplicate=True),
        Output('phenotype-data-store', 'data', allow_duplicate=True),
        Output('gwas-selected-label-store', 'data', allow_duplicate=True),
        Output('phenotype-label-store', 'data', allow_duplicate=True),
        #Output('pedigree-cytoscape', 'stylesheet', allow_duplicate=True),
    ],
    Input('selected-nodes-store', 'data'),
    [
        State('copy2gwas-selected-samples-store', 'data'),
        State('copy2pheno-selected-samples-store', 'data'),
        State('gwas-selected-samples-store', 'data'),
        State('phenotype-data-store', 'data'),
    ],
    prevent_initial_call=True
)
def sync_selected_nodes(selected_nodes, copy2_gwas_base, copy2_pheno_base,
                        prev_gwas_store, prev_pheno_store):
    from dash import no_update
    from dash.exceptions import PreventUpdate

    if not (selected_nodes or copy2_gwas_base or copy2_pheno_base):
        raise PreventUpdate

    nodes = selected_nodes or []
    base_gwas = copy2_gwas_base or []
    base_pheno = copy2_pheno_base or []

    # base phenotype id list
    base_pheno_ids = []
    for b in (base_pheno or []):
        if isinstance(b, dict) and "id" in b:
            base_pheno_ids.append(str(b["id"]))
        elif isinstance(b, str):
            base_pheno_ids.append(str(b))

    prev_gwas_store = prev_gwas_store or []
    if isinstance(prev_gwas_store, dict):
        prev_gwas_store = [prev_gwas_store]

    prev_pheno_store = prev_pheno_store or []
    if isinstance(prev_pheno_store, dict):
        prev_pheno_store = [prev_pheno_store]

    print(f"🔥 sync_selected_nodes 호출됨: 선택 {len(nodes)}개, GWAS기본 {len(base_gwas)}개, Pheno기본 {len(base_pheno)}개")

    # =====================================================
    # 🧬 1️⃣ GWAS 병합
    # =====================================================
    gwas_new = []
    gwas_label_dict = {}

    for nd in nodes:
        if nd.get("has_vcf") and nd.get("vcf_status"):
            vcf_id = str(nd["vcf_status"])
            name = (nd.get("id") or "").strip() or vcf_id
            gwas_new.append(vcf_id)
    
    base_gwas_ids = []
    base_gwas_dict = {}
    for item in base_gwas:
        if isinstance(item, str):
            vcf_id = item
            base_gwas_ids.append(vcf_id)
        elif isinstance(item, dict):
            vcf_id = str(item.get("vcf_status"))
            if not vcf_id:
                continue
            base_gwas_ids.append(vcf_id)
            if item.get("id"):
                base_gwas_dict[vcf_id] = item["id"]

    # base + new 병합
    gwas_seen = set()
    gwas_merged = []
    for val in base_gwas_ids  + gwas_new:
        if val and val not in gwas_seen:
            gwas_merged.append(val)
            gwas_seen.add(val)

    # 라벨 생성 (괄호 방식 + base 표시)
    for vcf_id in gwas_merged:
        # base 여부 표시
        is_base = vcf_id in (base_gwas or [])
        name = vcf_id
        for nd in nodes:
            if str(nd.get("vcf_status")) == vcf_id and nd.get("id"):
                name = nd["id"].strip()
                break
            else:
        # 2. base_gwas_dict에서 가져오기
                name = base_gwas_dict.get(vcf_id, vcf_id)
        label = f"{name} ({vcf_id}{', base' if is_base else ''})"
        gwas_label_dict[vcf_id] = label

    # 동일하면 no_update
    if prev_gwas_store == gwas_merged:
        print("⚖️ GWAS store 동일 → no_update 적용")
        gwas_store_output = no_update
        gwas_labels_output = no_update
    else:
        gwas_store_output = gwas_merged
        gwas_labels_output = gwas_label_dict

    # =====================================================
    # 🌾 2️⃣ Phenotype 병합
    # =====================================================
    pheno_new = []
    pheno_label_dict = {}

    for nd in nodes:
        if not nd.get("has_it"):
            continue

        it_number = None
        for key in ('it_number', 'id', 'variety_id', 'IT_Number'):
            val = nd.get(key)
            if val and str(val).strip() not in ('', 'N/A', 'None', 'null', 'No information'):
                it_number = str(val).strip()
                break

        if not it_number:
            continue

        name = (nd.get("id") or "").strip() or it_number
        pheno_new.append(it_number)

    # 병합
    pheno_seen = set()
    pheno_merged = []
    for val in base_pheno_ids + pheno_new:
        if val and val not in pheno_seen:
            pheno_merged.append(val)
            pheno_seen.add(val)

    pheno_store = [{'id': p} for p in pheno_merged]

    # 라벨 생성
    for pid in pheno_merged:
        is_base = pid in base_pheno_ids
        name = pid
        for nd in nodes:
            for key in ('it_number', 'id', 'variety_id', 'IT_Number'):
                if str(nd.get(key)) == pid and nd.get("id"):
                    name = nd["id"].strip()
                    break
        label = f"{name} ({pid}{', base' if is_base else ''})"
        pheno_label_dict[pid] = label

    if prev_pheno_store == pheno_store:
        print("⚖️ Phenotype store 동일 → no_update 적용")
        pheno_store_output = no_update
        pheno_labels_output = no_update
    else:
        print("⚖️ Phenotype store 변경 → update 적용")
        pheno_store_output = pheno_store
        pheno_labels_output = pheno_label_dict

    # =====================================================
    # 반환
    # =====================================================
    print(f"✅ GWAS {len(gwas_merged)}개 / PHENO {len(pheno_merged)}개 동기화 완료")
    base_stylesheet = get_default_stylesheet2()
    highlight_select_styles = []
    for nd in selected_nodes or []:
        nid = nd.get("id") or nd.get("name")
        if not nid:
            continue
        highlight_select_styles.append({
            "selector": f'node[id = "{nid}"]',
            "style": {
                "border-color": "#ff33cc",   # 짙은 초록 외곽선
                "border-width": 6,           # 두꺼운 강조선
                "transition-property": "border-color, border-width",
                "transition-duration": "0.3s",
            },
        })
    new_style = base_stylesheet + highlight_select_styles if highlight_select_styles else []

    #print(gwas_labels_output)
    return (
        gwas_store_output,
        pheno_store_output,
        gwas_labels_output,
        pheno_labels_output,
        #new_style
        
    )


@callback(
    Output('tab-content', 'children', allow_duplicate=True),
    [
        Input('analysis-tabs', 'active_tab'),
        State('fixed-vcf-item', 'data'),
        State('fixed-vcf-item-nopedi', 'data'),
    ],
    prevent_initial_call=True
)
def update_tab_view_gwas(active_tab, fixed_vcf_value, fixed_vcf_nopedi):
    """
    GWAS 탭 전환 시에만 반응.
    Phenotype 탭은 초기 로딩 시 이미 표시되어 있음.
    """
    from dash import no_update, html

    if active_tab != "gwas-tab":
        raise PreventUpdate

    print("🧬 GWAS 탭 전환 감지")

    vcf_status = None
    is_nopedi_case = isinstance(fixed_vcf_nopedi, dict) and fixed_vcf_nopedi.get("variety_id")
    is_regular_case = isinstance(fixed_vcf_value, dict) and fixed_vcf_value.get("processed_name")

    # status 추출
    if is_nopedi_case:
        vcf_status = fixed_vcf_nopedi.get("status")
    elif is_regular_case:
        vcf_status = fixed_vcf_value.get("status")

    print(f"vcf_status={vcf_status}")

    # GWAS 콘텐츠 생성
    if vcf_status and vcf_status != "No VCF":
        return create_gwas_analysis_content(vcf_status)

    # VCF 없음
    msg   = "No VCF data found for this nopedi variety." if is_nopedi_case else "VCF data is required for GWAS analysis."
    color = "warning" if is_nopedi_case else "danger"
    return html.Div([
        dbc.Alert([
            html.H5("❌ No VCF data", className="alert-heading"),
            html.P(msg)
        ], color=color)
    ])

@callback(
    Output('tab-content', 'children', allow_duplicate=True),
    [
        Input('analysis-tabs', 'active_tab'),
        Input('phenotype-data-store', 'data'),   # ✅ 사용자 선택 반응
    ],
    State('copy2pheno-selected-samples-store', 'data'),  # ✅ base 세트
    prevent_initial_call=True
)
def update_tab_view_phenotype(active_tab, pheno_store_data, copy2_pheno_base):
    """
    Phenotype 탭 전환 시 또는 phenotype-data-store 변경 시 패널 업데이트.
    - 탭 전환 시: base 세트(copy2pheno-selected-samples-store)로 렌더링
    - 데이터 변경 시: 선택 세트(phenotype-data-store)로 렌더링
    """
    from dash import html
    from dash.exceptions import PreventUpdate

    if active_tab != "phenotype-tab":
        raise PreventUpdate

    print("📊 Phenotype 탭 업데이트 감지")
    print(pheno_store_data)
    print(copy2_pheno_base)

    # ------------------------------
    # 1️⃣ 우선순위: phenotype-data-store → base 세트
    # ------------------------------
    if pheno_store_data:
        source = "store"
        raw_data = pheno_store_data
    else:
        source = "base"
        raw_data = copy2_pheno_base

    if not raw_data:
        return html.Div([
            dbc.Alert([
                html.H5("🌾 Phenotype Analysis Ready", className="alert-heading"),
                html.P("No phenotype data available.")
            ], color="info")
        ])

    # ------------------------------
    # 2️⃣ ID 추출
    # ------------------------------
    analysis_nodes = []
    for s in raw_data:
        if isinstance(s, dict):
            sid = s.get("id") or s.get("variety_id") or s.get("name")
            if sid:
                analysis_nodes.append(str(sid))
        elif isinstance(s, str):
            analysis_nodes.append(str(s))

    analysis_nodes = [x for x in analysis_nodes if x]
    analysis_nodes = list(dict.fromkeys(analysis_nodes))
    print(f"✅ Phenotype ({source}) 세트 노드: {analysis_nodes}")

    # ------------------------------
    # 3️⃣ PHENO_DF 유효성 검사
    # ------------------------------
    try:
        _empty = (PHENO_DF is None or PHENO_DF.empty)
    except NameError:
        _empty = True

    if _empty:
        return html.Div([
            dbc.Alert([
                html.H5("❌ No Phenotype Data", className="alert-heading"),
                html.P("PHENO_DF is empty or not loaded.")
            ], color="danger")
        ])

    # ------------------------------
    # 4️⃣ 패널 생성
    # ------------------------------
    return create_phenotype_panel_unified(
        selected_nodes=analysis_nodes,
        pheno_df=PHENO_DF,
        table_enabled=True
    )

# ─────────────────────────────────────
# 4️⃣ TOOLTIP 닫기: 주요 버튼 클릭 시 (reset/remove/expand/wide)
# ─────────────────────────────────────





@callback(
    Output("copy2selected-traits-store", "data",allow_duplicate=True),
    Input("integrated-selected-traits-store", "data"),
    State("group-type-radio", "value"),
    prevent_initial_call=True
)
def copy_selected_traits(integrated_traits, selected_type):
    """
    group-type이 'trait'일 때만 integrated-selected-traits-store → copy2selected-traits-store 로 복사
    """
    if selected_type != "trait":
        raise dash.exceptions.PreventUpdate

    integrated_traits = integrated_traits or []
    return integrated_traits


# ✅ 기존 “trait” 알람
make_alert_callbacks(
    store_id='integrated-selected-traits-store',
    span_id='new-trait-alert',
    timer_id='new-trait-timer',
    label_text='✨ New trait!',
    color='#e67e22'
)

# ✅ SNP presence (variants)
make_alert_callbacks(
    store_id='snp-occurrence-store',
    span_id='new-variant-alert',
    timer_id='new-variant-timer',
    label_text='🧬 New variant!',
    color='#16a085'
)

# ✅ GWAS samples
make_alert_callbacks(
    store_id='gwas-selected-samples-store',
    span_id='new-sample-alert',
    timer_id='new-sample-timer',
    label_text='🧫 New sample!',
    color='#2980b9'
)

@callback(
    Output('gwas-sample-scatter-wrapper', 'style'),
    Output('gwas-sample-table-wrapper', 'style'),
    Input('gwas-filter-tabs', 'value'),
    prevent_initial_call=True
)
def toggle_gwas_filter_tabs(tab_value):
    if tab_value == 'scatter':
        return {'display': 'block'}, {'display': 'none'}
    elif tab_value == 'table':
        return {'display': 'none'}, {'display': 'block'}
    return {'display': 'none'}, {'display': 'none'}

@callback(
    Output('gwas-sample-scatter-group-wrapper', 'style'),
    Output('gwas-sample-table-group-wrapper', 'style'),
    Input('gwas-group-tabs', 'value'),
    prevent_initial_call=True
)
def toggle_gwas_group_tabs(tab_value):
    if tab_value == 'Gscatter':
        return {'display': 'block'}, {'display': 'none'}
    elif tab_value == 'Gtable':
        return {'display': 'none'}, {'display': 'block'}
    return {'display': 'none'}, {'display': 'none'}

@callback(
    [
        Output('gwas-selected-badge', 'style', allow_duplicate=True),
        Output('gwas-mini-table-container', 'style', allow_duplicate=True),
        Output('gwas-mini-table-container', 'children', allow_duplicate=True),
        Output('gwas-selected-badge', 'children', allow_duplicate=True),
        Output('gwas-click-store', 'data', allow_duplicate=True),  # 유지 (초기화용)
    ],
    [
        Input('gwas-selected-badge', 'n_clicks'),
        Input('selected-nodes-store', 'data'),
        Input('fixed-vcf-item', 'data'),
        Input('fixed-vcf-item-nopedi', 'data'),
        Input('gwas-click-store', 'data'),  # ✅ click_store를 Input으로 승격
    ],
     [
        State('gwas-selected-label-store', 'data'),  # ✅ 추가
    ],
    prevent_initial_call=True
)
def toggle_badge_and_table(
    n_clicks, selected_nodes, fixed_vcf_item, fixed_vcf_item_nopedi, click_store,gwas_label_dict  
):
    from dash import html, ctx, no_update
    from dash.exceptions import PreventUpdate
    import re
    gwas_label_dict = gwas_label_dict or {}
    #print(gwas_label_dict)
    #print(click_store)
    # ---------------------------------------------
    # 내부 유틸 함수
    # ---------------------------------------------
    def remove_parentheses_numbers(st):
        if not st:
            return st
        return re.sub(r'\s*\(\d+\)$', '', str(st)).strip()

    def get_subtrait_color(subtrait, color_map=None):
        """Subtrait에 대한 색상 반환"""
        if not color_map:
            ALL_SUBTRAITS = [
                'Yield', 'Stress', 'Plantvigor', 'Biochemical',
                'Plantgrowth_Development', 'Plant_quality',
                'Biologicalprocess', 'Plantmorphology', 'Sterility_Fertility'
            ]
            try:
                import plotly.express as px
                ALL_SUBTRAIT_COLORS = px.colors.qualitative.Set1[:9]
                color_map = {sub: color for sub, color in zip(ALL_SUBTRAITS, ALL_SUBTRAIT_COLORS)}
                color_map['default'] = '#888888'
            except ImportError:
                color_map = {'default': '#888888'}
        return color_map.get(subtrait, color_map.get('default', '#888888'))

    trigger = ctx.triggered_id
    click_store = click_store or []
    click_count = len(click_store)
    badge_label = f"Selected Variants ({click_count})"

    # ---------------------------------------------
    # 1️⃣ node 변경 → badge 리셋 & container 닫기
    # ---------------------------------------------
    if trigger in ['selected-nodes-store', 'fixed-vcf-item', 'fixed-vcf-item-nopedi']:
        print("🔄 Node changed → reset badge & hide container")
        return (
            {
                "border": "none",
                "padding": "2px 6px",
                "cursor": "pointer",
                "transition": "all 0.2s ease",
            },
            {"display": "none"},
            no_update,
            "Selected Variants (0)",
            [],
        )

    # ---------------------------------------------
    # 2️⃣ badge 클릭 이벤트
    # ---------------------------------------------
    if trigger == 'gwas-selected-badge':
        # badge 토글 시 열기/닫기 제어
        if not click_store:
            return (
                {
                    "border": "none",
                    "padding": "2px 6px",
                    "cursor": "pointer",
                    "transition": "all 0.2s ease",
                },
                {"display": "none"},
                no_update,
                badge_label,
                click_store,
            )

        # 홀수 클릭 → 열기
        if n_clicks % 2 == 1:
            badge_style = {
                "border": "2px solid #f39c12",
                "borderRadius": "8px",
                "padding": "2px 6px",
                "cursor": "pointer",
                "transition": "all 0.2s ease",
            }

            spans = []
            for item in click_store:
                vid = item.get("Variant_ID", "")
                chr_ = item.get("Chromosome", "")
                pos = item.get("Position", "")
                trait = item.get("Trait", "")
                samples = item.get("Samples", [])
                subtrait = item.get("Subtrait", "")

                if not vid:
                    continue

                # ✅ Position 변환
                try:
                    pos_val = int(float(pos))
                    pos_display = f"{pos_val/1_000_000:.2f}M" if pos_val >= 1_000_000 else f"{pos_val:,}"
                except:
                    pos_display = str(pos)

                # ✅ GT 표시 (3개 이상은 ...)
                #gt_pairs = [f"{s}: {item.get(f'{s}_GT', '')}" for s in samples if item.get(f"{s}_GT")]
                #gt_display = ", ".join(gt_pairs[:3]) + ", ..." if len(gt_pairs) > 3 else ", ".join(gt_pairs)
                gt_pairs = [
                    f"{gwas_label_dict.get(s, s)}: {item.get(f'{s}_GT', '')}"
                    for s in samples if item.get(f"{s}_GT")
                ]
                gt_display = ", ".join(gt_pairs[:3]) + ", ..." if len(gt_pairs) > 3 else ", ".join(gt_pairs)
                tooltip_text = f"chr{chr_}:{pos_display}\nTrait: {trait}\n" + "\n".join(gt_pairs)
                short_text = f"chr{chr_}:{pos_display} | {trait} | {gt_display}"

                # ✅ 색상
                subtrait_clean = remove_parentheses_numbers(subtrait)
                subtrait_color = get_subtrait_color(subtrait_clean)

                spans.append(
                    html.Span(
                        short_text,
                        id={"type": "variant-span", "index": vid},
                        n_clicks=0,
                        title=tooltip_text,
                        className="variant-span",
                        style={
                            "margin": "4px 5px",
                            "padding": "4px 10px",
                            "border": f"1.5px solid {subtrait_color}",
                            "borderRadius": "8px",
                            "display": "inline-block",
                            "cursor": "pointer",
                            "backgroundColor": f"{subtrait_color}20",
                            "color": 'black',
                            "fontSize": "12px",
                            "fontWeight": "600",
                            "boxShadow": "1px 1px 3px rgba(0,0,0,0.1)",
                            "userSelect": "none",
                            "transition": "all 0.2s ease",
                        },
                    )
                )

            container = html.Div(
                spans,
                style={
                    "display": "flex",
                    "flexWrap": "wrap",
                    "gap": "6px",
                    "padding": "6px 0",
                    "backgroundColor": "#f8f9fa",
                    "border": "1px solid #dee2e6",
                    "borderRadius": "8px",
                    "paddingLeft": "10px",
                    "paddingRight": "10px",
                    "boxShadow": "inset 0 1px 2px rgba(0,0,0,0.05)",
                },
            )
            return badge_style, {"display": "block"}, container, badge_label, click_store

        # 짝수 클릭 → 닫기
        else:
            return (
                {
                    "border": "none",
                    "padding": "2px 6px",
                    "cursor": "pointer",
                    "transition": "all 0.2s ease",
                },
                {"display": "none"},
                no_update,
                badge_label,
                click_store,
            )

    # ---------------------------------------------
    # 3️⃣ click_store 업데이트 시 자동 갱신
    # ---------------------------------------------
    if trigger == 'gwas-click-store':
        print("♻️ click_store updated → refreshing badge (if open)")

        # badge가 닫혀 있다면 container 갱신 안 함
        # 여기서 display 상태는 기억되어 있지 않으므로, 항상 열기 기준으로 재구성 가능
        if not click_store:
            return (
                {
                    "border": "none",
                    "padding": "2px 6px",
                    "cursor": "pointer",
                    "transition": "all 0.2s ease",
                },
                {"display": "none"},
                no_update,
                "Selected Variants (0)",
                click_store,
            )

        # 열린 상태일 때 자동 refresh
        badge_style = {
            "border": "2px solid #f39c12",
            "borderRadius": "8px",
            "padding": "2px 6px",
            "cursor": "pointer",
            "transition": "all 0.2s ease",
        }

        spans = []
        for item in click_store:
            vid = item.get("Variant_ID", "")
            chr_ = item.get("Chromosome", "")
            pos = item.get("Position", "")
            trait = item.get("Trait", "")
            samples = item.get("Samples", [])
            subtrait = item.get("Subtrait", "")

            if not vid:
                continue

            try:
                pos_val = int(float(pos))
                pos_display = f"{pos_val/1_000_000:.2f}M" if pos_val >= 1_000_000 else f"{pos_val:,}"
            except:
                pos_display = str(pos)
            
            gt_pairs = [
                    f"{gwas_label_dict.get(s, s)}: {item.get(f'{s}_GT', '')}"
                    for s in samples if item.get(f"{s}_GT")
                ]

            #gt_pairs = [f"{s}: {item.get(f'{s}_GT', '')}" for s in samples if item.get(f"{s}_GT")]
            gt_display = ", ".join(gt_pairs[:3]) + ", ..." if len(gt_pairs) > 3 else ", ".join(gt_pairs)

            subtrait_clean = remove_parentheses_numbers(subtrait)
            subtrait_color = get_subtrait_color(subtrait_clean)

            short_text = f"chr{chr_}:{pos_display} | {trait} | {gt_display}"

            spans.append(
                html.Span(
                    short_text,
                    id={"type": "variant-span", "index": vid},
                    n_clicks=0,
                    className="variant-span",
                    style={
                        "margin": "4px 5px",
                        "padding": "4px 10px",
                        "border": f"1.5px solid {subtrait_color}",
                        "borderRadius": "8px",
                        "display": "inline-block",
                        "cursor": "pointer",
                        "backgroundColor": f"{subtrait_color}20",
                        "color": 'black',
                        "fontSize": "12px",
                        "fontWeight": "600",
                        "boxShadow": "1px 1px 3px rgba(0,0,0,0.1)",
                        "userSelect": "none",
                        "transition": "all 0.2s ease",
                    },
                )
            )

        container = html.Div(spans, style={
            "display": "flex",
            "flexWrap": "wrap",
            "gap": "6px",
            "padding": "6px 0",
            "backgroundColor": "#f8f9fa",
            "border": "1px solid #dee2e6",
            "borderRadius": "8px",
            "paddingLeft": "10px",
            "paddingRight": "10px",
            "boxShadow": "inset 0 1px 2px rgba(0,0,0,0.05)",
        })

        return badge_style, {"display": "block"}, container, badge_label, click_store

    raise PreventUpdate

@callback(
    [
        Output('pedigree-cytoscape', 'stylesheet', allow_duplicate=True),
        Output({'type': 'variant-span', 'index': ALL}, 'children', allow_duplicate=True),
        Output({'type': 'variant-span', 'index': ALL}, 'style', allow_duplicate=True),
        Output('gwas-genotype-legend', 'children', allow_duplicate=True),
    ],
    [Input({'type': 'variant-span', 'index': ALL}, 'n_clicks')],
    [
        State('gwas-click-store', 'data'),
        State('selected-nodes-store', 'data'),
        State('fixed-vcf-item', 'data'),
        State('fixed-vcf-item-nopedi', 'data'),
        State('highlighted-components-store', 'data'),
        State('component-info-store', 'data'),
        State('gwas-selected-label-store', 'data'),
    ],
    prevent_initial_call=True
)
def highlight_from_variant_spans(
    n_clicks_list, click_store, selected_nodes,
    fixed_vcf_item, fixed_vcf_item_nopedi,
    selected_comps, comp_info, gwas_label_dict
):
    import re
    from dash import ctx, no_update, html
    from dash.exceptions import PreventUpdate

    gwas_label_dict = gwas_label_dict or {}
    click_store = click_store or []
    selected_nodes = selected_nodes or []
    if not n_clicks_list:
        raise PreventUpdate

    triggered = ctx.triggered_id
    if not triggered or "index" not in triggered:
        raise PreventUpdate
    variant_id = triggered["index"]

    # -------------------------------
    # Helper functions
    # -------------------------------
    def contains_minor(gt, minor):
        if not gt or not minor:
            return False
        gt = str(gt).upper().replace("|", "/").strip()
        minor = str(minor).upper().strip()
        alleles = [a.strip() for a in re.split(r"[\/|]", gt) if a.strip()]
        return minor in alleles

    def remove_parentheses_numbers(st):
        return re.sub(r"\s*\(\d+\)$", "", str(st)).strip() if st else st

    def get_subtrait_color(subtrait):
        import plotly.express as px
        ALL_SUBTRAITS = [
            'Yield', 'Stress', 'Plantvigor', 'Biochemical',
            'Plantgrowth_Development', 'Plant_quality',
            'Biologicalprocess', 'Plantmorphology', 'Sterility_Fertility'
        ]
        color_map = {s: c for s, c in zip(ALL_SUBTRAITS, px.colors.qualitative.Set1[:9])}
        color_map['default'] = '#888888'
        return color_map.get(subtrait, color_map['default'])

    # -------------------------------
    # Merge selected + fixed nodes
    # -------------------------------
    fixed_nodes = []
    if fixed_vcf_item:
        fixed_nodes.append({
            "id": fixed_vcf_item.get("processed_name", ""),
            "vcf_id": fixed_vcf_item.get("status", "")
        })
    if fixed_vcf_item_nopedi:
        fixed_nodes.append({
            "id": fixed_vcf_item_nopedi.get("variety_id", ""),
            "vcf_id": fixed_vcf_item_nopedi.get("status", "")
        })
    all_nodes = selected_nodes + fixed_nodes

    # -------------------------------
    # Determine mode (single/multi)
    # -------------------------------
    click_counts = n_clicks_list or []
    variant_idx = next((i for i, v in enumerate(click_store)
                        if v.get("Variant_ID") == variant_id), None)
    if variant_idx is None:
        raise PreventUpdate

    click_count = click_counts[variant_idx] if len(click_counts) > variant_idx else 0
    clicked_active = (click_count % 2 == 1)
    any_other_active = any((c or 0) % 2 == 1 for i, c in enumerate(click_counts) if i != variant_idx)

    if clicked_active:
        mode = "single"
    elif not clicked_active and not any_other_active:
        mode = "multi"
    elif not clicked_active and any_other_active:
        mode = "single"
    else:
        mode = "multi"

    # -------------------------------
    # Determine target variants
    # -------------------------------
    if mode == "single":
        target_variants = [v for v in click_store if v["Variant_ID"] == variant_id]
    else:
        target_variants = click_store[:] if click_store else []

    # -------------------------------
    # Compute node colors
    # -------------------------------
    node_colors = {}
    for node in all_nodes:
        nid = str(node.get("id", ""))
        vcf_id = str(node.get("vcf_id", "")).strip()
        if not vcf_id:
            continue

        matched_flags = []
        for var in target_variants:
            samples = var.get("Samples", [])
            sample_match = next((s for s in samples if s == vcf_id or s in vcf_id or vcf_id in s), None)
            if not sample_match:
                continue
            gt = var.get(f"{sample_match}_GT", "")
            minor = var.get("Minor_Allele", "")
            matched_flags.append(contains_minor(gt, minor))

        if not matched_flags:
            node_colors[nid] = "#323232"
        elif all(matched_flags):
            node_colors[nid] = "#9B111E"
        elif any(matched_flags):
            node_colors[nid] = "#f1c40f"
        else:
            node_colors[nid] = "#323232"


    # -------------------------------
    # Update span children & styles
    # -------------------------------
    span_children, span_styles = [], []
    for item, n_click in zip(click_store, n_clicks_list):
        vid = item.get("Variant_ID")
        if not vid:
            span_children.append(no_update)
            span_styles.append(no_update)
            continue

        subtrait_clean = remove_parentheses_numbers(item.get("Subtrait", item.get("Trait", "")))
        subtrait_color = get_subtrait_color(subtrait_clean)
        samples = item.get("Samples", [])
        gt_pairs = [
            f"{gwas_label_dict.get(s, s)}: {item.get(f'{s}_GT', '')}"
            for s in samples if item.get(f"{s}_GT")
        ]
        chr_, pos, trait = item.get("Chromosome", ""), item.get("Position", ""), item.get("Trait", "")

        try:
            pos_val = int(float(pos))
            pos_display = f"{pos_val/1_000_000:.2f}M" if pos_val >= 1_000_000 else f"{pos_val:,}"
        except:
            pos_display = str(pos)

        is_clicked = (vid == variant_id and clicked_active)

        # ✅ GT 표시 규칙
        if len(gt_pairs) <= 3:
            gt_text = ", ".join(gt_pairs)  # 3개 이하 → 항상 전체 표시
        elif is_clicked:
            gt_text = ", ".join(gt_pairs)  # 클릭된 항목 → 전체 표시
        else:
            gt_text = ", ".join(gt_pairs[:3]) + ", ..."  # 기본 → 앞 3개 + ellipsis

        # ✅ 텍스트 구성
        text = f"chr{chr_}:{pos_display} | {trait} | {gt_text}"

        # ✅ 스타일 정의
        style_dict = {
            "margin": "4px 5px",
            "padding": "4px 10px",
            "border": f"1.5px solid {subtrait_color}",
            "borderRadius": "8px",
            "display": "inline-block",
            "cursor": "pointer",
            "backgroundColor": f"{subtrait_color}50" if is_clicked else f"{subtrait_color}20",
            "color": 'black',
            "fontSize": "12px",
            "fontWeight": "600",
            "boxShadow": (
                "2px 2px 5px rgba(0,0,0,0.25)" if is_clicked
                else "1px 1px 3px rgba(0,0,0,0.1)"
            ),
            "userSelect": "none",
            "transition": "all 0.2s ease",
            "whiteSpace": "nowrap",
            "overflow": "hidden",
            "textOverflow": "clip" if is_clicked or len(gt_pairs) <= 3 else "ellipsis",
            "maxWidth": "none" if is_clicked or len(gt_pairs) <= 3 else "600px",
        }

        span_children.append(text)
        span_styles.append(style_dict)

    # -------------------------------
    # Cytoscape style
    # -------------------------------
    base_styles = get_default_stylesheet2()
    overlay_styles = [{
        "selector": f'[id = "{nid}"]',
        "style": {"background-color": color, "border-color": "#f39c12", "border-width": "3px"},
    } for nid, color in node_colors.items()]
    new_styles = base_styles + overlay_styles

    # -------------------------------
    # Legend
    # -------------------------------
    legend_children = html.Div([
        html.Div(["● ", "All minor alleles present"], style={"color": "#9B111E", "fontSize": "12px"}),
        html.Div(["● ", "Partially matched"], style={"color": "#f1c40f", "fontSize": "12px"}),
        html.Div(["● ", "No minor allele"], style={"color": "#323232", "fontSize": "12px"}),
        html.Div(f"Mode: {'Single variant' if mode == 'single' else 'Multi variant'}",
                 style={"fontSize": "12px", "fontStyle": "italic", "color": "#555"})
    ])

    return new_styles, span_children, span_styles, legend_children







'''
@callback(
    Output('pedigree-cytoscape', 'stylesheet', allow_duplicate=True),
    Input('gwas-mini-table', 'active_cell'),
    State('gwas-mini-table', 'data'),
    State('gwas-click-store', 'data'),
    State('gwas-selected-samples-store', 'data'),
    State('available-nodes-store', 'data'),
    State('addv-apply-groups-nopedi', 'data'),
    prevent_initial_call=True
)
def highlight_from_badge_table(active_cell, table_data, click_store,
                               selected_samples, avail_nodes, group_nopedi):
    """Badge mini-table에서 클릭 시 pedigree 강조"""
    import re
    from dash.exceptions import PreventUpdate

    if not active_cell or not table_data:
        raise PreventUpdate

    row_idx = active_cell.get("row")
    if not (0 <= row_idx < len(table_data)):
        raise PreventUpdate

    row = table_data[row_idx]
    variant_id = str(row.get("Variant_ID"))
    if not variant_id:
        raise PreventUpdate

    print(f"🎯 [BADGE→PEDIGREE] variant={variant_id}")

    click_store = click_store or []
    record = next((v for v in click_store if str(v["Variant_ID"]) == variant_id), None)
    if not record:
        raise PreventUpdate

    minor = str(record.get("Minor_Allele") or "")
    selected_samples = selected_samples or []
    avail_nodes = avail_nodes or []
    nodes_np = (group_nopedi or {}).get("records", []) if isinstance(group_nopedi, dict) else []

    sample_gt_map = {s: record.get(f"{s}_GT", "") for s in selected_samples}

    node_map = {}
    for n in (avail_nodes + nodes_np):
        if n.get("vcf_id") and n.get("id"):
            node_map[str(n["vcf_id"])] = str(n["id"])

    node_gt_map = {
        node_map.get(s): gt
        for s, gt in sample_gt_map.items()
        if s in node_map
    }

    node_colors = {}
    for nid, gt in node_gt_map.items():
        if not gt:
            node_colors[nid] = "#bcbcbc"
        elif minor and re.search(rf"\b{re.escape(minor)}\b", gt, re.IGNORECASE):
            node_colors[nid] = "#27ae60"
        else:
            node_colors[nid] = "#7f8c8d"

    base_styles = get_default_stylesheet2()
    overlay_styles = [
        {
            "selector": f'node[id="{nid}"]',
            "style": {
                "background-color": color,
                "border-color": color,
                "border-width": "3px",
            },
        }
        for nid, color in node_colors.items() if nid
    ]

    return base_styles + overlay_styles
    '''
'''
@callback(
    Output('scatter-refresh-trigger', 'data', allow_duplicate=True),
    Input('gwas-click-store', 'data'),
    prevent_initial_call=True
)
def trigger_scatter_refresh(click_store):
    """gwas-click-store가 비면 scatter 재실행 트리거"""
    if not click_store or len(click_store) == 0:
        import time
        return time.time()  # ✅ timestamp 사용 → 항상 값 변경
    raise PreventUpdate
'''

@callback(
    [
        Output('profiles-panel', 'style'),
        Output('profiles-content', 'children')
    ],
    Input('btn-profiles-toggle', 'n_clicks'),
    State('available-nodes-store', 'data'),
    prevent_initial_call=False
)
def toggle_profiles_panel(n_clicks, available_nodes):
    """
    📋 Detailed Profiles 패널 토글 + 내용 표시
    - non-pedigree는 별도 처리하지 않음 (사각형 노드로 이미 포함됨)
    """
    from dash import no_update, html
    import dash_bootstrap_components as dbc

    nodes = available_nodes or []

    # --- 패널 표시 제어
    show_panel = bool(n_clicks and n_clicks % 2 == 1)
    panel_style = {
        "display": "block" if show_panel else "none",
        "backgroundColor": "#f8f9fa",
        "border": "1px solid #dee2e6",
        "borderRadius": "6px",
        "marginTop": "8px",
        "padding": "10px",
        "maxHeight": "220px",
        "overflowY": "auto",
        "boxShadow": "inset 0 1px 3px rgba(0,0,0,0.1)",
    }

    # --- 내용 렌더링
    if not nodes:
        content = html.Div([
            html.Div("📋 No profiles to display", style={
                "color": "#6c757d",
                "fontStyle": "italic",
                "textAlign": "center",
                "fontSize": "13px"
            })
        ])
    else:
        content = []
        for node in nodes:
            name = node.get("name") or node.get("id", "Unknown")
            itnum = node.get("it_number", "N/A")
            vcf = "✅" if node.get("has_vcf") else "❌"
            pheno = "✅" if node.get("has_pheno") else "❌"
            shape = "rectangle" if node.get("type") == "nopedi" else "circle"

            content.append(html.Div([
                html.H6(name, style={
                    "margin": "0 0 4px 0",
                    "color": "#2c3e50",
                    "fontWeight": "bold"
                }),
                html.Div([
                    html.Span("IT: ", style={"fontWeight": "bold"}),
                    html.Span(itnum), html.Br(),
                    html.Span("VCF: ", style={"fontWeight": "bold"}),
                    html.Span(vcf), html.Br(),
                    html.Span("Phenotype: ", style={"fontWeight": "bold"}),
                    html.Span(pheno), html.Br(),
                    html.Span("Shape: ", style={"fontWeight": "bold"}),
                    html.Span("Non-pedigree (rectangle)" if shape == "rectangle" else "Pedigree (circle)"),
                ], style={"fontSize": "12px", "lineHeight": "1.4"})
            ], style={
                "backgroundColor": "white",
                "padding": "6px 8px",
                "marginBottom": "8px",
                "borderRadius": "6px",
                "border": "1px solid #e9ecef",
                "boxShadow": "0 1px 2px rgba(0,0,0,0.05)"
            }))

    return panel_style, content



clientside_callback(
    """
    function(overData){
      const tip = document.getElementById('tooltip');
      if(!tip){ return null; }

      if(!overData){
        return null; // hide는 assets의 스크립트가 처리
      }

      window.__lastHoverTick = Date.now();

      const val = (o, k, d) => (o && o[k] !== undefined && o[k] !== null ? o[k] : d);
      const asArr = (v) => Array.isArray(v) ? v : (v==null || v==='' ? [] : [String(v)]);
      const esc = (s) => String(s).replace(/[&<>"]/g, m => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[m]));

      const variety_id    = esc(val(overData, 'id', 'Unknown'));
      const variety_label = esc(val(overData, 'label', variety_id));

      const parents = asArr(val(overData, 'parents', []));
      const children= asArr(val(overData, 'children', []));
      const it_num  = esc(val(overData, 'it_number', 'No information'));
      const vcf_stat= esc(val(overData, 'vcf_status', 'No VCF'));
      const has_it  = !!val(overData, 'has_it', false);
      const has_vcf = !!val(overData, 'has_vcf', false);

      function summarize(list, singular, plural){
        if (!list.length) return `No ${plural}`;
        const head = list.slice(0,2).map(esc);
        if (list.length > 2) return `${list.length} ${plural}: ${head.join(', ')}...`;
        return `${list.length} ${list.length===1?singular:plural}: ${head.join(', ')}`;
      }
      const parents_info  = summarize(parents, 'parent', 'parents');
      const children_info = summarize(children, 'child', 'children');
      const phenotype_info= (has_it && it_num!=='No information') ? `✅ IT: ${it_num}` : '❌ No IT number';
      const vcf_info      = (has_vcf && vcf_stat!=='No VCF')       ? `✅ VCF: ${vcf_stat}` : '❌ No VCF data';

      const html = `
        <div style="font-size:12px">
          <!-- Header -->
          <div style="display:flex;align-items:flex-start;border-bottom:1px solid #555;padding-bottom:8px;margin-bottom:8px">
            <div style="flex:1">
              <strong style="color:#fff;font-size:16px;display:block;margin-bottom:6px">${variety_label}</strong>
              <small style="color:#ccc;font-size:11px">ID: ${variety_id}</small>
            </div>
            <span
              title="Tooltip"
              style="font-size:18px;opacity:0.7;line-height:1;user-select:none">×</span>
          </div>

          <!-- Pedigree info -->
          <div style="margin-bottom:8px">
            <div style="margin-bottom:4px">
              <i class="fas fa-level-up-alt" style="margin-right:6px;color:#ffc107"></i>
              <span style="font-size:11px">${esc(parents_info)}</span>
            </div>
            <div style="margin-bottom:4px">
              <i class="fas fa-level-down-alt" style="margin-right:6px;color:#17a2b8"></i>
              <span style="font-size:11px">${esc(children_info)}</span>
            </div>
          </div>

          <!-- Data availability -->
          <div>
            <div style="margin-bottom:4px">
              <i class="fas fa-chart-bar" style="margin-right:6px;color:#28a745"></i>
              <span style="font-size:11px">${esc(phenotype_info)}</span>
            </div>
            <div style="margin-bottom:6px">
              <i class="fas fa-dna" style="margin-right:6px;color:#6f42c1"></i>
              <span style="font-size:11px">${esc(vcf_info)}</span>
            </div>
          </div>

          <!-- Footer -->
          <div>
            <i class="fas fa-mouse-pointer" style="margin-right:4px;color:#6c757d"></i>
            <span style="font-size:9px;font-style:italic;color:#aaa">Click to expand/select</span>
          </div>
        </div>
      `;

      tip.innerHTML   = html;
      tip.style.display    = 'block';
      tip.style.visibility = 'visible';
      tip.style.opacity    = '1';
      tip.style.transform  = 'translateZ(0)';

      return 1; // 더미 Store 갱신
    }
    """,
    Output('tooltip-sync', 'data'),
    Input('pedigree-cytoscape', 'mouseoverNodeData'),
)


@callback(
    Output('gwas-selected-label-store', 'data',allow_duplicate=True),
    Input('copy2gwas-selected-samples-store', 'data'),
    prevent_initial_call='initial_duplicate'  # ✅ 앱 시작 시에도 실행되도록
)
def initialize_gwas_label_store(base_gwas):
    gwas_label_dict = {}

    base_gwas = base_gwas or []
    for item in base_gwas:
        if isinstance(item, str):
            vcf_id = item
            name = vcf_id
        elif isinstance(item, dict):
            vcf_id = str(item.get("vcf_status"))
            name = item.get("id", vcf_id)
        else:
            continue

        label = f"{name} ({vcf_id}, base)"
        gwas_label_dict[vcf_id] = label

    return gwas_label_dict



@callback(
    [
        Output("trait-legend-section", "style", allow_duplicate=True),
        Output("gwas-genotype-legend-section", "style", allow_duplicate=True),
    ],
    Input("analysis-tabs", "active_tab"),
    prevent_initial_call=True
)
def toggle_legend_sections(active_tab):
    base_style = {
        "padding": "15px",
        "backgroundColor": "#ffffff",
        "borderRadius": "8px",
        "border": "1px solid #dee2e6",
        "marginBottom": "15px",
    }

    if active_tab == "gwas-tab":
        return (
            {"display": "none", **base_style},
            {"display": "block", **base_style},
        )
    else:
        return (
            {"display": "block", **base_style},
            {"display": "none", **base_style},
        )


@callback(
    Output('group-type-radio', 'options'),
    [
        Input('integrated-selected-traits-store', 'data'),
        Input('presence-sets-store', 'data'),
        Input('gwas-selected-samples-store', 'data'),
    ],
)
def update_group_type_labels(traits_data, presence_data, samples_data):
    """
    각 store 값의 길이에 따라 count를 표시하되,
    0이면 count를 숨기고 기본 라벨만 유지
    """
    # --- count 계산 ---
    traits_count = len(traits_data or [])
    presence_count = len(presence_data or []) if isinstance(presence_data, list) else (1 if presence_data else 0)
    samples_count = len(samples_data or [])

    # --- count를 조건부 표시 ---
    def fmt(label, count):
        return f"{label} ({count})" if count > 0 else label

    options = [
        {'label': fmt("Trait", traits_count), 'value': 'trait'},
        {'label': fmt("Variant ID", presence_count), 'value': 'variant'},
        {'label': fmt("Sample", samples_count), 'value': 'sample'},
    ]
    return options