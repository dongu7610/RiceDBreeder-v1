# app_csv_backend.py
import dash
from dash import html, dcc, Input, Output, State, callback, ALL, MATCH, ctx
import dash_bootstrap_components as dbc
import plotly.graph_objs as go
import pandas as pd
import numpy as np
import dash_table
from datetime import datetime
from assets.pedi_utils import create_header, create_footer, COMMON_STYLES

# ===== CSV 백엔드 함수들 =====
# data_backend_csv.py가 같은 폴더에 있어야 합니다.
from data_backend_csv import (
    get_db_data,               # resource_types에 따라 phenotype + resource join (CSV)
    get_available_columns,     # 카테고리 라벨 구성 (하드매핑)
    get_column_type_from_db,   # trait 타입 판정 (fallback)
    get_resource_types,        # 리소스 타입 목록 (CSV)
    build_results_table_data   # 대형 SQL 대체(VCF/thal3 매칭 포함)
)

# ===== App =====
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

# Pre-allocate filter IDs
TOTAL_FILTERS = 61
FILTER_IDS = [f"filter_{i}" for i in range(1, TOTAL_FILTERS + 1)]

# ==== FilterManager (현 구조 유지) ====
class FilterManager:
    def __init__(self):
        self.used_filters = set()
        self.filter_sequence = []
        self.current_index = 0

    def get_next_available_id(self):
        """Returns both the filter_id string and numeric index"""
        self.current_index += 1
        filter_id = f"filter_{self.current_index}"
        self.used_filters.add(filter_id)
        self.filter_sequence.append(filter_id)
        return filter_id, self.current_index

    def remove_filter(self, filter_id):
        """Remove a filter and update the sequence"""
        if filter_id in self.used_filters:
            self.used_filters.remove(filter_id)
            self.filter_sequence = [f for f in self.filter_sequence if f != filter_id]
            self.reindex_filters()

    def reindex_filters(self):
        new_sequence = []
        new_used_filters = set()
        for i, _old in enumerate(self.filter_sequence, 1):
            new_id = f"filter_{i}"
            new_sequence.append(new_id)
            new_used_filters.add(new_id)
        self.filter_sequence = new_sequence
        self.used_filters = new_used_filters
        self.current_index = len(self.filter_sequence)

    def get_filter_sequence(self):
        return self.filter_sequence.copy()

filter_manager = FilterManager()

# ===== Helper functions (DB 의존 제거판) =====
def get_initial_filter_values(df, column_type, column):
    try:
        if column_type == 'date':
            numeric_values = pd.to_numeric(df[column], errors='coerce').dropna()
            if not numeric_values.empty:
                min_date = int(numeric_values.min())
                max_date = int(numeric_values.max())
                # mmdd 해석 안전 처리
                def _mm(x, dflt): 
                    s = str(int(x))
                    return int(s[:-2]) if len(s) >= 3 else dflt
                def _dd(x, dflt):
                    s = str(int(x))
                    return int(s[-2:]) if len(s) >= 2 else dflt
                return {
                    'date_range': [min_date, max_date],
                    'start_month': _mm(min_date, 1),
                    'start_day': _dd(min_date, 1),
                    'end_month': _mm(max_date, 12),
                    'end_day': _dd(max_date, 31)
                }
        elif column_type == 'numeric':
            numeric_values = pd.to_numeric(df[column], errors='coerce').dropna()
            if not numeric_values.empty:
                return {'range': [float(numeric_values.min()), float(numeric_values.max())]}
        elif column_type == 'categorical':
            unique_values = df[column].dropna().unique()
            if len(unique_values) > 0:
                return {'categories': list(unique_values)}
    except Exception:
        pass

    if column_type == 'date':
        return {'date_range': [101, 1231], 'start_month': 1, 'start_day': 1, 'end_month': 12, 'end_day': 31}
    if column_type == 'numeric':
        return {'range': [0, 100]}
    if column_type == 'categorical':
        return {'categories': []}
    return {}

def create_filter_visualization(original_df, input_df, filtered_df, column, column_type):
    try:
        fig = go.Figure()
        base_count = len(filtered_df)
        input_count = len(input_df)

        if column_type == 'date':
            all_values = pd.to_numeric(input_df[column], errors='coerce').dropna()
            filtered_values = pd.to_numeric(filtered_df[column], errors='coerce').dropna()
            if len(all_values) == 0:
                return html.Div("No valid date values found in the data")

            months = []
            for x in all_values:
                try:
                    s = str(int(x))
                    if len(s) >= 3:
                        m = int(s[:-2])
                        if 1 <= m <= 12:
                            months.append(m)
                except Exception:
                    continue
            months = sorted(set(months))
            if not months:
                return html.Div("Could not extract valid months from the date values")

            date_ranges, date_labels, x_positions = [], [], []
            for month in months:
                ranges = [
                    (float(f"{month}01"), float(f"{month}10")),
                    (float(f"{month}11"), float(f"{month}20")),
                    (float(f"{month}21"), float(f"{month}31"))
                ]
                date_ranges.extend(ranges)
                date_labels.extend([f"{month} 1-10", f"{month} 11-20", f"{month} 21-last day"])

            all_counts, filtered_counts = [], []
            for i, (start, end) in enumerate(date_ranges):
                all_counts.append(len(all_values[(all_values >= start) & (all_values <= end)]))
                filtered_counts.append(len(filtered_values[(filtered_values >= start) & (filtered_values <= end)]))
                x_positions.append(i)

            fig.add_trace(go.Bar(x=x_positions, y=all_counts, name='data', opacity=0.5))
            fig.add_trace(go.Bar(x=x_positions, y=filtered_counts, name='filterd data', opacity=0.7))
            fig.update_layout(
                xaxis=dict(ticktext=date_labels, tickvals=x_positions, tickangle=45, title="interval"),
                yaxis=dict(title="Counts"),
                barmode='overlay',
                showlegend=True
            )

        elif column_type == 'numeric':
            numeric_values = pd.to_numeric(input_df[column], errors='coerce')
            filtered_values = pd.to_numeric(filtered_df[column], errors='coerce')
            fig.add_trace(go.Histogram(x=numeric_values, name="data", opacity=0.5))
            fig.add_trace(go.Histogram(x=filtered_values, name="filterd data", opacity=0.7))
            fig.update_layout(barmode='overlay', showlegend=True)

        else:  # categorical
            value_counts = input_df[column].value_counts()
            filtered_counts = filtered_df[column].value_counts()
            fig.add_trace(go.Bar(x=value_counts.index.astype(str), y=value_counts.values, name='data', opacity=0.5))
            if not filtered_df.empty:
                fig.add_trace(go.Bar(x=filtered_counts.index.astype(str), y=filtered_counts.values, name='filterd data', opacity=0.7))

        fig.update_layout(title=f"{column}", height=300, margin={'t': 30, 'b': 30, 'l': 30, 'r': 30}, showlegend=True)

        return html.Div([
            html.Strong(
                f"Filtered results: {base_count:,} out of total {input_count:,} cases",
                className="filter-result-text",
                style={'fontSize': '14px','color': '#2c3e50','display': 'block','marginBottom': '10px',
                       'marginTop': '15px','padding': '8px','backgroundColor': '#f8f9fa',
                       'borderRadius': '4px','border': '1px solid #dee2e6'}
            ),
            dcc.Graph(figure=fig),
        ])
    except Exception as e:
        return f"Error creating visualization: {str(e)}"

def get_filter_key(column, column_type):
    return f"{column}_{column_type}"

def create_filter_component(filter_id, column, column_type, df):
    filter_label = f"Filter : {column}"
    months = [{'label': f'{i}', 'value': i} for i in range(1, 13)]
    days = [{'label': f'{i}', 'value': i} for i in range(1, 32)]
    initial_values = get_initial_filter_values(df, column_type, column)

    # numeric range
    if column_type == 'numeric':
        range_vals = initial_values.get('range', [0, 100])
        range_min = float(range_vals[0]); range_max = float(range_vals[1])
        range_step = (range_max - range_min) / 100 if range_min != range_max else 1
    else:
        range_min = 0; range_max = 100; range_step = 1

    if column_type == 'date':
        initial_start_month = initial_values.get('start_month')
        initial_start_day = initial_values.get('start_day')
        initial_end_month = initial_values.get('end_month')
        initial_end_day = initial_values.get('end_day')
    else:
        initial_start_month = initial_start_day = initial_end_month = initial_end_day = None

    return html.Div([
        html.Div([
            html.H4(filter_label, style={'display':'inline-block','fontSize':'16px','color':'#2c3e50','fontWeight':'bold','margin':'0'}),
            html.Button('×', id={'type': 'remove-filter', 'index': filter_id},
                        style={'float':'right','border':'none','background':'none','fontSize':'20px','cursor':'pointer','color':'#dc3545','fontWeight':'bold'})
        ], style={'marginBottom': '15px'}),

        # Date
        html.Div([
            html.Div([
                html.Label('Start date:', style={'fontWeight':'bold','marginRight':'20px','width':'80px','display':'inline-block','verticalAlign':'top','paddingTop':'5px','fontSize':'14px','color':'#495057'}),
                html.Div([
                    html.Div([html.Label('Month:'), dcc.Dropdown(id={'type':'filter-startdate-month', 'index': filter_id}, options=months, value=initial_start_month, placeholder='month', style={'width':'100px'})],
                             style={'display': 'inline-block','marginRight': '10px'}),
                    html.Div([html.Label('Day:'), dcc.Dropdown(id={'type':'filter-startdate-day', 'index': filter_id}, options=days, value=initial_start_day, placeholder='day', style={'width':'100px'})],
                             style={'display': 'inline-block'})
                ])
            ], style={'marginBottom': '10px'}),

            html.Div([
                html.Label('End date:', style={'fontWeight':'bold','marginRight':'20px','width':'80px','display':'inline-block','verticalAlign':'top','paddingTop':'5px','fontSize':'14px','color':'#495057'}),
                html.Div([
                    html.Div([html.Label('Month:'), dcc.Dropdown(id={'type':'filter-enddate-month', 'index': filter_id}, options=months, value=initial_end_month, placeholder='month', style={'width':'100px'})],
                             style={'display': 'inline-block','marginRight': '10px'}),
                    html.Div([html.Label('Day:'), dcc.Dropdown(id={'type':'filter-enddate-day', 'index': filter_id}, options=days, value=initial_end_day, placeholder='day', style={'width':'100px'})],
                             style={'display': 'inline-block'})
                ])
            ])
        ], id={'type': 'date-container', 'index': filter_id}, style={'display': 'block' if column_type == 'date' else 'none'}),

        # Numeric
        html.Div([
            dcc.RangeSlider(
                id={'type': 'numeric-range', 'index': filter_id},
                min=range_min, max=range_max, step=range_step, value=[range_min, range_max],
                tooltip={'placement': 'bottom', 'always_visible': True},
                marks={range_min: f'{range_min:.2f}', range_max: f'{range_max:.2f}'}
            )
        ], id={'type': 'numeric-container', 'index': filter_id}, style={'display': 'block' if column_type == 'numeric' else 'none'}),

        # Categorical
        html.Div([
            dcc.Dropdown(
                id={'type': 'categorical-filter', 'index': filter_id},
                options=[{'label': str(val), 'value': str(val)} for val in sorted(df[column].dropna().unique()) if pd.notna(val)],
                multi=True, placeholder="Select categories...",
                value=[str(val) for val in sorted(df[column].dropna().unique()) if pd.notna(val)]
            )
        ], id={'type': 'categorical-container', 'index': filter_id}, style={'display': 'none' if column_type != 'categorical' else 'block'}),

        html.Div(id={'type': 'filter-result', 'index': filter_id}),
        html.Hr()
    ], className="filter-block")

def process_id_list(df):
    if df.empty:
        return []
    return sorted(df['id'].unique().tolist())

# ===== Layout =====
layout = html.Div([
    # Store components
    dcc.Store(id='filter-store', data={'filters': []}),
    dcc.Store(id='available-columns-store', data=[]),
    dcc.Store(id='filter-chain-store', data={'filters': [], 'results': []}),
    dcc.Store(id='filter-update-trigger', data=0),
    dcc.Store(id='last-filter-result', data=None),  # 마지막 필터 결과 저장용

    create_header(),

    html.Div([
        # Left sidebar with filters
        html.Div([
            html.Div([
                html.H2("Filtering Section",
                        style={'marginBottom':'10px','color':'#2c3e50','borderBottom':'2px solid #3498db','paddingBottom':'8px'}),
                html.Div(["Search the entire rice resource database.", html.Br(),
                          "Select a Resource Type and click Add Filter to customize your search."],
                         style={'fontSize':'14px','color':'#666','marginBottom':'20px','lineHeight':'1.5'})
            ], className="filter-section-title"),

            html.Div([
                # Resource Type Selection
                html.Div([
                    html.Div([
                        html.Div("Resource Type Selection", className="filter-title"),
                        dcc.Dropdown(
                            id='filter-x',
                            options=[{'label': rt, 'value': rt} for rt in get_resource_types()],
                            placeholder="Select Resource Types",
                            multi=True
                        ),
                        html.Div(id='filter-x-result', className="filter-content"),
                    ], className="filter-content"),
                ], className="filter-block"),

                # Add Filter Button
                html.Div([
                    dbc.Button("Add filter(+)", id="add-filter-btn", color="primary", size="sm",
                               className="add-filter-btn", style={'display': 'none'}),
                ], className="filter-btn-container"),

                # Dynamic Filters Container
                html.Div(id='dynamic-filters', className="scrollable-filters"),
            ], className="filter-container"),
        ], className="sidebar", style={'height': '100%', 'borderRight': '1px solid #dee2e6'}),

        # Right content area
        html.Div([
            html.Div([
                html.Div([
                    html.H2("Table Section",
                            style={'marginBottom':'10px','color':'#2c3e50','borderBottom':'2px solid #3498db','paddingBottom':'8px'}),
                    html.Div(["View multiple filtering options in the table and navigate to a specific variety card via its link."],
                             style={'fontSize':'14px','color':'#666','marginBottom':'20px','lineHeight':'1.5'})
                ], className="result-section-title"),
                html.Div(id='result-table-container', className="result-content")
            ], className="result-section"),
        ], className="main-content"),
    ], className="content-container"),

    # Filter Selection Modal
    dbc.Modal([
        dbc.ModalHeader("Choose Filter"),
        dbc.ModalBody([
            dcc.Dropdown(id='filter-type-selector', placeholder="Select a filter to add")
        ]),
        dbc.ModalFooter([
            dbc.Button("add", id="confirm-filter-btn", className="ms-auto", color="primary"),
            dbc.Button("cancel", id="close-modal-btn", className="ms-2")
        ]),
    ], id="filter-modal"),

    create_footer()
], style={
    'backgroundColor': COMMON_STYLES['background_color'],
    'fontFamily': COMMON_STYLES['font_family'],
    'minHeight': '100vh',
    'display': 'flex',
    'flexDirection': 'column'
})

app.layout = layout

# ===== Callbacks =====
@callback(
    Output('filter-chain-store', 'data'),
    [Input({'type': 'filter-startdate-month', 'index': ALL}, 'value'),
     Input({'type': 'filter-startdate-day', 'index': ALL}, 'value'),
     Input({'type': 'filter-enddate-month', 'index': ALL}, 'value'),
     Input({'type': 'filter-enddate-day', 'index': ALL}, 'value'),
     Input({'type': 'numeric-range', 'index': ALL}, 'value'),
     Input({'type': 'categorical-filter', 'index': ALL}, 'value')],
    [State('filter-chain-store', 'data')]
)
def update_filter_chain(start_months, start_days, end_months, end_days,
                        range_values, categorical_values, current_chain):
    if not current_chain:
        current_chain = {'filters': [], 'results': []}
    ctx_triggered = ctx.triggered[0] if ctx.triggered else None
    if not ctx_triggered:
        return current_chain

    try:
        trigger_id = ctx_triggered['prop_id'].split('.')[0]
        if not trigger_id:
            return current_chain

        trigger_info = eval(trigger_id)
        filter_type = trigger_info['type']
        filter_index = trigger_info['index']

        filter_to_update = next((f for f in current_chain['filters'] if f['id'] == filter_index), None)
        if not filter_to_update:
            return current_chain

        if filter_type.startswith('filter-startdate') or filter_type.startswith('filter-enddate'):
            if all([start_months[filter_index-1], start_days[filter_index-1],
                    end_months[filter_index-1], end_days[filter_index-1]]):
                start_date = start_months[filter_index-1] * 100 + start_days[filter_index-1]
                end_date = end_months[filter_index-1] * 100 + end_days[filter_index-1]
                filter_to_update['filter_values'] = {
                    'date_range': [start_date, end_date],
                    'start_month': start_months[filter_index-1],
                    'start_day': start_days[filter_index-1],
                    'end_month': end_months[filter_index-1],
                    'end_day': end_days[filter_index-1]
                }

        elif filter_type == 'numeric-range':
            if range_values[filter_index-1]:
                filter_to_update['filter_values'] = {'range': range_values[filter_index-1]}

        elif filter_type == 'categorical-filter':
            if categorical_values[filter_index-1]:
                filter_to_update['filter_values'] = {'categories': categorical_values[filter_index-1]}

        return current_chain
    except Exception:
        return current_chain

@callback(
    Output('filter-x-result', 'children'),
    [Input('filter-x', 'value')]
)
def update_filter_x_result(selected_resources):
    if not selected_resources:
        return "No resource types selected"
    try:
        df = get_db_data(selected_resources)
        if df.empty:
            return "No data found for selected types"
        id_list = process_id_list(df)
        return html.Div([html.Div(f"Total IDs: {len(id_list)}", className="result-stats")])
    except Exception as e:
        return f"Error loading data: {str(e)}"

@callback(
    Output('add-filter-btn', 'style'),
    [Input('filter-x', 'value')]
)
def toggle_filter_button(selected_resources):
    if selected_resources and len(selected_resources) > 0:
        return {'display': 'block'}
    return {'display': 'none'}

@callback(
    [Output('dynamic-filters', 'children', allow_duplicate=True),
     Output('filter-chain-store', 'data', allow_duplicate=True),
     Output('filter-modal', 'is_open', allow_duplicate=True),
     Output('filter-type-selector', 'options', allow_duplicate=True)],
    [Input('add-filter-btn', 'n_clicks'),
     Input('confirm-filter-btn', 'n_clicks'),
     Input('close-modal-btn', 'n_clicks'),
     Input({'type': 'remove-filter', 'index': ALL}, 'n_clicks'),
     Input('filter-x', 'value')],
    [State('filter-type-selector', 'value'),
     State('dynamic-filters', 'children'),
     State('filter-chain-store', 'data')],
    prevent_initial_call=True
)
def manage_filters(add_clicks, confirm_clicks, close_clicks, remove_clicks,
                   resource_type_value, selected_column, existing_filters, filter_chain):
    ctx_triggered = ctx.triggered_id
    if existing_filters is None:
        existing_filters = []
    if filter_chain is None:
        filter_chain = {'filters': [], 'results': []}

    try:
        # 사용 가능한 컬럼 옵션
        if resource_type_value:
            df = get_db_data(resource_type_value)
            all_columns = get_available_columns(df)
            used_columns = {f['column'] for f in filter_chain.get('filters', [])}
            available_options = [col for col in all_columns if col['value'] not in used_columns]
        else:
            available_options = []

        # 필터 제거
        if isinstance(ctx_triggered, dict) and ctx_triggered.get('type') == 'remove-filter':
            target_filter_id = ctx_triggered['index']
            filter_manager.remove_filter(f"filter_{target_filter_id}")

            updated_filters, new_chain_filters = [], []
            for i, filter_component in enumerate(existing_filters, 1):
                if isinstance(filter_component, dict) and 'props' in filter_component:
                    header_div = next(
                        (child for child in filter_component['props']['children']
                         if isinstance(child, dict) and 'props' in child and 'children' in child['props']
                         and isinstance(child['props']['children'], list) and len(child['props']['children']) > 1
                         and isinstance(child['props']['children'][1], dict)
                         and child['props']['children'][1].get('props', {}).get('id', {}).get('index') == target_filter_id),
                        None
                    )
                    if not header_div:
                        updated_component = update_filter_component_ids(filter_component, i)
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

        if ctx_triggered == 'add-filter-btn':
            return existing_filters, filter_chain, True, available_options

        if ctx_triggered == 'close-modal-btn':
            return existing_filters, filter_chain, False, available_options

        if ctx_triggered == 'confirm-filter-btn' and selected_column:
            df = get_db_data(resource_type_value)
            next_filter_id, filter_number = filter_manager.get_next_available_id()
            if next_filter_id is None:
                raise Exception("Maximum number of filters reached")

            column_type = get_column_type_from_db(selected_column)
            new_filter = create_filter_component(filter_number, selected_column, column_type, df)
            initial_filter_values = get_initial_filter_values(df, column_type, selected_column)

            filter_chain['filters'].append({
                'id': filter_number,
                'column': selected_column,
                'type': column_type,
                'sequence': len(filter_manager.filter_sequence),
                'filter_values': initial_filter_values
            })
            return existing_filters + [new_filter], filter_chain, False, available_options

        return existing_filters, filter_chain, False, available_options
    except Exception:
        return existing_filters, filter_chain, False, available_options

def update_filter_component_ids(component, new_index):
    if not isinstance(component, dict) or 'props' not in component:
        return component
    updated_component = component.copy()
    props = updated_component['props']
    if 'children' in props:
        props['children'] = [update_child_component_ids(child, new_index) for child in props['children']]
    return updated_component

def update_child_component_ids(child, new_index):
    if not isinstance(child, dict) or 'props' not in child:
        return child
    updated_child = child.copy()
    props = updated_child['props']
    if 'id' in props and isinstance(props['id'], dict):
        props['id'] = {'type': props['id']['type'], 'index': new_index}
    if 'children' in props:
        if isinstance(props['children'], list):
            props['children'] = [update_child_component_ids(g, new_index) for g in props['children']]
        elif isinstance(props['children'], dict):
            props['children'] = update_child_component_ids(props['children'], new_index)
    return updated_child

@callback(
    Output({'type': 'filter-result', 'index': MATCH}, 'children'),
    [Input({'type': 'filter-startdate-month', 'index': MATCH}, 'value'),
     Input({'type': 'filter-startdate-day', 'index': MATCH}, 'value'),
     Input({'type': 'filter-enddate-month', 'index': MATCH}, 'value'),
     Input({'type': 'filter-enddate-day', 'index': MATCH}, 'value'),
     Input({'type': 'numeric-range', 'index': MATCH}, 'value'),
     Input({'type': 'categorical-filter', 'index': MATCH}, 'value'),
     Input('filter-x', 'value'),
     Input('filter-update-trigger', 'data')],
    [State({'type': 'filter-result', 'index': MATCH}, 'id'),
     State('filter-chain-store', 'data')]
)
def update_filter_results(start_month, start_day, end_month, end_day,
                          range_value, categorical_value,
                          resource_type_value, _trigger,
                          result_id, filter_chain):
    if not resource_type_value:
        return "No data selected"
    try:
        original_df = get_db_data(resource_type_value)
        current_df = original_df.copy()

        current_filter_id = result_id['index']
        if isinstance(current_filter_id, str) and '_' in current_filter_id:
            current_filter_id = int(current_filter_id.split('_')[1])

        current_filter = next((f for f in filter_chain['filters'] if f['id'] == current_filter_id), None)
        if not current_filter:
            return "Filter not found in chain"

        all_filters = sorted(filter_chain['filters'], key=lambda x: x['sequence'])
        current_sequence = current_filter['sequence']

        for f in all_filters:
            if f['sequence'] >= current_sequence:
                break
            current_df = apply_stored_filter(current_df, f)

        input_df = current_df.copy()
        column = current_filter['column']
        column_type = current_filter['type']

        if column_type == 'date':
            numeric_values = pd.to_numeric(current_df[column], errors='coerce')
            current_df = current_df[numeric_values.notna()]
            if all([start_month, start_day, end_month, end_day]):
                start_date = start_month * 100 + start_day
                end_date = end_month * 100 + end_day
                numeric_values = pd.to_numeric(current_df[column], errors='coerce')
                mask = (numeric_values >= start_date) & (numeric_values <= end_date)
                current_df = current_df[mask]

        elif column_type == 'numeric':
            numeric_values = pd.to_numeric(current_df[column], errors='coerce')
            current_df = current_df[numeric_values.notna()]
            if range_value:
                current_df = current_df[numeric_values.between(*range_value)]

        elif column_type == 'categorical':
            current_df = current_df[current_df[column].notna()]
            if categorical_value:
                current_df = current_df[current_df[column].isin(categorical_value)]

        return create_filter_visualization(original_df, input_df, current_df, column, column_type)
    except Exception as e:
        return f"Error applying filter: {str(e)}"

def apply_stored_filter(df, filter_info):
    filtered_df = df.copy()
    column = filter_info['column']
    ftype = filter_info['type']
    fvals = filter_info.get('filter_values', {})

    try:
        if ftype == 'date' and 'date_range' in fvals:
            s, e = fvals['date_range']
            vals = pd.to_numeric(filtered_df[column], errors='coerce')
            mask = vals.notna() & (vals >= s) & (vals <= e)
            filtered_df = filtered_df[mask]
        elif ftype == 'numeric' and 'range' in fvals:
            rmin, rmax = fvals['range']
            vals = pd.to_numeric(filtered_df[column], errors='coerce')
            filtered_df = filtered_df[vals.between(rmin, rmax)]
        elif ftype == 'categorical' and 'categories' in fvals:
            filtered_df = filtered_df[filtered_df[column].isin(fvals['categories'])]
    except Exception:
        return df
    return filtered_df

@callback(
    Output('result-table-container', 'children'),
    [Input({'type': 'filter-result', 'index': ALL}, 'children'),
     Input('filter-chain-store', 'data'),
     Input('add-filter-btn', 'n_clicks')],
    [State('filter-x', 'value')]
)
def update_results_table(_filter_results, filter_chain, _add_clicks, resource_type_value):
    if not resource_type_value:
        return "Please select resource types"
    if not filter_chain or 'filters' not in filter_chain or not filter_chain['filters']:
        return "Add filters to see results"

    try:
        df = get_db_data(resource_type_value)

        def format_filter_condition(filter_info):
            column = filter_info['column']
            ftype = filter_info['type']
            fvals = filter_info.get('filter_values', {})
            if ftype == 'date':
                if 'date_range' in fvals:
                    s, e = fvals['date_range']
                    return f"{column}: {s}-{e}"
                return f"{column}: date filter"
            elif ftype == 'numeric':
                if 'range' in fvals:
                    rmin, rmax = fvals['range']
                    return f"{column}: {rmin:.1f}-{rmax:.1f}"
                return f"{column}: numeric filter"
            elif ftype == 'categorical':
                if 'categories' in fvals:
                    cats = fvals['categories']
                    return f"{column}: {', '.join(cats[:3])}" + (f"... (+{len(cats)-3} more)" if len(cats) > 3 else "")
                return f"{column}: categorical filter"
            return f"{column}: unknown filter"

        all_filters = sorted(filter_chain['filters'], key=lambda x: x['sequence'])
        for f in all_filters:
            df = apply_stored_filter(df, f)

        filtered_ids = df['id'].tolist() if not df.empty else []
        if not filtered_ids:
            individual = [f"Filter{i}: {format_filter_condition(fi)}" for i, fi in enumerate(all_filters, 1)]
            return html.Div([
                html.H4("No data found after applying filters", className="table-title"),
                html.Div([
                    html.Strong("Applied Filters:", style={'color':'#2c3e50','marginBottom':'8px','display':'block'}),
                    html.Div([html.Div(cond, style={'fontSize':'14px','color':'#6c757d','marginBottom':'4px','paddingLeft':'10px'}) for cond in individual])
                ], style={'textAlign':'left','margin':'20px 0','padding':'15px','backgroundColor':'#f8f9fa','borderRadius':'5px','border':'1px solid #dee2e6'})
            ])

        # === 여기서부터 대형 SQL 대체 ===
        result_df = build_results_table_data(filtered_ids)  # CSV 기반 합성(VCF/thal3 포함)
        # 각 필터 컬럼의 '실제 값' 추가(원본 df에서 끌어오기)
        original_df = get_db_data(resource_type_value)

        filter_columns = []
        for i, f in enumerate(all_filters, 1):
            col = f['column']
            out_col = f"Filter{i}_{col}"
            filter_columns.append(out_col)
            id_to_val = dict(zip(original_df['id'], original_df[col]))
            result_df[out_col] = result_df['IT Number'].map(id_to_val)

        # 필터 조건 텍스트
        individual = [f"Filter{i}: {format_filter_condition(fi)}" for i, fi in enumerate(all_filters, 1)]

        base_columns = [
            {'name': 'IT Number', 'id': 'IT Number', 'type': 'text', 'presentation': 'markdown'},
            {'name': 'Rice Variety Name', 'id': 'Rice Variety Name'},
            {'name': 'Rice Line Name', 'id': 'Rice Line Name'},
            {'name': 'Variety_type', 'id': 'Variety_type'},
            {'name': 'VCF Status', 'id': 'VCF Status'},
            {'name': 'Pedigree Status', 'id': 'detection_status'}
        ]
        add_columns = [{'name': f"Filter{i+1}: {all_filters[i]['column']}", 'id': fc, 'type': 'text'} for i, fc in enumerate(filter_columns)]
        all_columns = base_columns + add_columns

        table = dash_table.DataTable(
            data=[{
                **row,
                'IT Number': f'[{row["IT Number"]}](/app7_2?id={row["IT Number"]}&processed_name={row["processed_name"]})'
            } for row in result_df.to_dict('records')],
            columns=all_columns,
            page_size=20,
            page_action='native',
            sort_action='native',
            sort_mode='multi',
            filter_action='native',
            markdown_options={'link_target': '_blank'},
            style_cell={'fontSize': '14px'},
            style_header={'fontSize': '14px', 'fontWeight': 'bold'},
            style_data_conditional=[
                {'if': {'column_id': c['id']}, 'fontSize': '12px', 'color': '#0066cc', 'textAlign': 'center', 'fontWeight': '500'}
                for c in add_columns
            ]
        )

        return html.Div([
            html.H4(f"Multi-Filtering Results: Total {len(result_df):,}", className="table-title"),
            html.Div([
                html.Strong("Applied Filters:", style={'color':'#2c3e50','marginBottom':'8px','display':'block'}),
                html.Div([html.Div(cond, style={'fontSize':'14px','color':'#6c757d','marginBottom':'4px','paddingLeft':'10px'}) for cond in individual])
            ], style={'textAlign':'left','margin':'10px 0 20px 0','padding':'15px','backgroundColor':'#f8f9fa','borderRadius':'5px','border':'1px solid #dee2e6'}),
            table
        ], className="table-container")
    except Exception as e:
        return f"Error displaying results: {str(e)}"

@callback(
    Output('filter-update-trigger', 'data'),
    [Input({'type': 'filter-startdate-month', 'index': ALL}, 'value'),
     Input({'type': 'filter-startdate-day', 'index': ALL}, 'value'),
     Input({'type': 'filter-enddate-month', 'index': ALL}, 'value'),
     Input({'type': 'filter-enddate-day', 'index': ALL}, 'value'),
     Input({'type': 'numeric-range', 'index': ALL}, 'value'),
     Input({'type': 'categorical-filter', 'index': ALL}, 'value')]
)
def update_filter_trigger(*_args):
    return datetime.now().timestamp()

@callback(
    Output('last-filter-result', 'data'),
    [Input({'type': 'filter-result', 'index': ALL}, 'children'),
     Input('filter-chain-store', 'data')]
)
def update_last_filter_result(filter_results, filter_chain):
    if not filter_chain or 'filters' not in filter_chain or not filter_chain['filters']:
        return None
    last_filter = max(filter_chain['filters'], key=lambda x: x['sequence'])
    last_filter_id = last_filter['id']
    for result in filter_results:
        if isinstance(result, dict) and 'props' in result:
            filter_id = result['props'].get('id', {}).get('index')
            if filter_id == last_filter_id:
                return result
    return None

@callback(
    Output({'type': 'categorical-filter', 'index': MATCH}, 'value'),
    [Input({'type': 'categorical-filter', 'index': MATCH}, 'value')],
    [State('filter-x', 'value'),
     State('filter-chain-store', 'data'),
     State({'type': 'categorical-filter', 'index': MATCH}, 'id')]
)
def handle_categorical_selection(selected_values, resource_type_value, filter_chain, component_id):
    if not ctx.triggered:
        filter_index = component_id['index']
    else:
        trigger = ctx.triggered[0]
        if not trigger:
            return selected_values
        trigger_id = eval(trigger['prop_id'].split('.')[0])
        filter_index = trigger_id['index']

    try:
        current_filter = next((f for f in filter_chain['filters'] if f['id'] == filter_index), None)
        if (not selected_values or len(selected_values) == 0) and current_filter and resource_type_value:
            if current_filter['type'] == 'categorical':
                df = get_db_data(resource_type_value)
                column = current_filter['column']
                all_values = sorted([str(val) for val in df[column].dropna().unique() if pd.notna(val)])
                return all_values
    except Exception:
        pass
    return selected_values

@callback(
    Output('dynamic-filters', 'style'),
    Input('filter-chain-store', 'data')
)
def adjust_filter_container_height(filter_chain):
    num_filters = len(filter_chain.get('filters', [])) if filter_chain and 'filters' in filter_chain else 0
    base_height = 400
    additional_height_per_filter = 300
    max_additional_filters = 3
    if num_filters <= max_additional_filters:
        dynamic_height = base_height + (num_filters * additional_height_per_filter)
        return {'minHeight': f'{dynamic_height}px', 'maxHeight': f'{dynamic_height + 100}px', 'overflowY': 'visible'}
    else:
        max_height = base_height + (max_additional_filters * additional_height_per_filter)
        return {'minHeight': f'{max_height}px', 'maxHeight': f'{max_height}px', 'overflowY': 'auto'}

if __name__ == '__main__':
    app.run_server(debug=True, port=8082)
