# app_home.py  — CSV backend only (no DB calls)
from dash import html, dcc, callback
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import dash
import pandas as pd

# 기존 공용 UI 유틸은 그대로 사용
from assets.pedi_utils import create_header, create_footer, COMMON_STYLES

# 👉 홈 전용 CSV 헬퍼
from backend_csv_home import build_variety_df

# ------------------------------
# 스타일
# ------------------------------
CARD_STYLE = {
    'height': '375px',
    'cursor': 'pointer',
    'transition': 'transform 0.2s',
    'marginBottom': '20px',
    'border': '2px solid #dee2e6',
}

# ------------------------------
# 데이터 준비 (CSV → 단 1회 로드)
# ------------------------------
VARIETY_DF: pd.DataFrame = build_variety_df()

def _to_options(df: pd.DataFrame):
    if df is None or df.empty:
        return []
    return [
        {
            'label': row['display_name'],
            'value': row['id'],
            'vcf_status': row['vcf_status'],
            'has_pedigree': bool(row['has_pedigree']),
            'pedigree_name': (row['pedigree_name'] if pd.notna(row['pedigree_name']) else None),
        }
        for _, row in df.iterrows()
    ]

INITIAL_OPTIONS = _to_options(VARIETY_DF)

# ------------------------------
# 공용 UI 조각
# ------------------------------
def create_info_div(row: pd.Series):
    """품종 정보 DIV 생성"""
    return html.Div([
        html.Div([
            html.Div([
                html.Span("VCF Status: ", style={'fontWeight': 'bold'}),
                html.Span(row['vcf_status'],
                          style={'color': 'green' if row['vcf_status'] != 'No VCF' else 'red'})
            ], style={'marginBottom': '5px'}),
            html.Div([
                html.Span("Pedigree Status: ", style={'fontWeight': 'bold'}),
                html.Span("Available" if row['has_pedigree'] else "Not Available",
                          style={'color': 'green' if row['has_pedigree'] else 'red'})
            ])
        ], style={
            'padding': '10px',
            'backgroundColor': '#f8f9fa',
            'borderRadius': '5px',
            'fontSize': '1.1rem'
        })
    ])

# ------------------------------
# 레이아웃
# ------------------------------
layout = html.Div([
    create_header(),

    # 🔹 Header 아래 Tutorial 버튼 (우측 상단)
    html.Div([
        html.Div([
            dcc.Link([
                html.I(className="fas fa-book-open", style={
                    'fontSize': '28px',
                    'color': '#4A90E2',
                    'marginRight': '6px'
                }),
                html.Span("Tutorial", style={
                    'color': '#4A90E2',
                    'fontWeight': '600',
                    'fontSize': '16px',
                    'verticalAlign': 'middle'
                })
            ],
            href='/tutorial',
            target='_blank',
            style={
                'textDecoration': 'none',
                'display': 'flex',
                'alignItems': 'center',
                'justifyContent': 'flex-end',
                'gap': '6px'
            })
        ],
        style={
            'display': 'flex',
            'justifyContent': 'flex-end',
            'margin': '10px 60px 0 0'
        })
    ]),

    html.Div([
        dbc.Container([
            dbc.Row([
                # 필터 페이지로 이동 카드 (app7_1은 그대로 유지)
                dbc.Col([

                    dbc.Card([
                        dbc.CardBody([
                            html.H3("Phenotype Search", className="card-title mb-4 text-center"),

                            html.P("Search accessions by phenotype filters.",
                                className="text-muted text-center mb-4",
                                style={'fontSize': '1.1rem'}),

                            html.Div([
                                dcc.Link([
                                    html.Div([
                                        html.I(className="fas fa-filter", style={'fontSize': '48px', 'color': '#0066cc'}),
                                        html.P("Go to filtering page", style={'marginTop': '10px', 'color': '#0066cc', 'fontWeight': 'bold'})
                                    ], className="text-center")
                                ],
                                href='/app7_1',
                                target="_blank",
                                style={'textDecoration': 'none'})
                            ])
                        ], className="p-5 text-center")
                    ], style=CARD_STYLE, className="h-100 shadow-sm hover-lift")

                ], width=12, md=6),

                # 품종 검색 카드 (CSV 기반)
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H3("Rice Variety Search", className="card-title mb-4"),
                            html.Div([
                                dcc.Dropdown(
                                    id='itnumber-input',
                                    options=INITIAL_OPTIONS,
                                    placeholder='Enter or select an IT Number',
                                    style={'width': '100%', 'marginBottom': '20px', 'fontSize': '1.1rem'},
                                    searchable=True,
                                    clearable=True
                                )
                            ], style={'minHeight': '38px'}),

                            dcc.Loading(
                                id="loading-info",
                                type="default",
                                children=[
                                    html.Div(id='selected-variety-info',
                                             className="my-4",
                                             style={'minHeight': '20px'})
                                ]
                            ),

                            html.A(
                                dbc.Button(
                                    "Go to variety card",
                                    id='search-button',
                                    color="primary",
                                    className="w-100",
                                    style={'height': '50px', 'fontSize': '1.1rem'}
                                ),
                                id='variety-link',
                                href='',
                                target="_blank",
                                style={'textDecoration': 'none', 'width': '100%', 'display': 'block'}
                            )
                        ], className="p-5")
                    ], style=CARD_STYLE, className="h-100 shadow-sm")
                ], width=12, md=6),
            ], className="py-4"),
        ]),
    ], style={
        'flex': '1',
        'display': 'flex',
        'alignItems': 'center',
    }),

    create_footer()
], style={
    'backgroundColor': COMMON_STYLES['background_color'],
    'fontFamily': COMMON_STYLES['font_family'],
    'minHeight': '100vh',
    'display': 'flex',
    'flexDirection': 'column'
})

# ------------------------------
# 콜백
# ------------------------------
@callback(
    Output('variety-link', 'href'),
    [Input('itnumber-input', 'value')],
    [State('itnumber-input', 'options')]
)
def update_variety_link(search_value, options):
    if not search_value:
        return '#'
    # 선택된 옵션에서 정보 조회
    selected = next((opt for opt in (options or []) if opt['value'] == search_value), None)
    if not selected:
        return '#'
    processed_name = 'nopedi' if not selected.get('has_pedigree') else selected.get('pedigree_name')
    return f'/app7_2?id={search_value}&processed_name={processed_name}'

@callback(
    Output('itnumber-input', 'options'),
    Input('itnumber-input', 'search_value'),
    prevent_initial_call=True
)
def update_options(search_value):
    """드롭다운 search_value로 로컬 필터 (CSV 재조회 없음)"""
    if not search_value:
        return INITIAL_OPTIONS
    needle = str(search_value).lower()
    return [opt for opt in INITIAL_OPTIONS if needle in str(opt['label']).lower()]

@callback(
    [Output('selected-variety-info', 'children'),
     Output('selected-variety-info', 'style')],
    Input('itnumber-input', 'value'),
    prevent_initial_call=True
)
def update_variety_info(selected_value):
    base_style = {'minHeight': '60px', 'transition': 'all 0.3s ease'}

    if not selected_value:
        return "", {**base_style, 'visibility': 'hidden'}
    
    row = VARIETY_DF.loc[VARIETY_DF['id'] == selected_value]
    if row.empty:
        return "", {**base_style, 'visibility': 'hidden'}
    
    return create_info_div(row.iloc[0]), {**base_style, 'visibility': 'visible'}






#