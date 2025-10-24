#multi_filter.py
# multi_filter.py
# 멀티 필터 UI/로직 모듈화: 기존 함수들을 클래스화해서 재사용
from typing import Dict, List, Any
import pandas as pd
from dash import html, dcc
import plotly.graph_objects as go
from datetime import datetime
from data_backend_csv import determine_trait_type_fallback

class FilterManager:
    def __init__(self):
        self.used_filters = set()
        self.filter_sequence: List[str] = []
        self.current_index = 0

    def get_next_available_id(self):
        self.current_index += 1
        filter_id = f"filter_{self.current_index}"
        self.used_filters.add(filter_id)
        self.filter_sequence.append(filter_id)
        return filter_id, self.current_index

    def remove_filter(self, filter_id: str):
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

def get_filter_key(column: str, column_type: str) -> str:
    return f"{column}_{column_type}"

def get_initial_filter_values(df: pd.DataFrame, column_type: str, column: str) -> Dict[str, Any]:
    try:
        if column_type == 'date':
            numeric_values = pd.to_numeric(df[column], errors='coerce').dropna()
            if not numeric_values.empty:
                min_date = int(numeric_values.min())
                max_date = int(numeric_values.max())
                return {
                    'date_range': [min_date, max_date],
                    'start_month': int(str(min_date)[:-2]) if min_date >= 100 else 1,
                    'start_day': int(str(min_date)[-2:]) if min_date >= 100 else 1,
                    'end_month': int(str(max_date)[:-2]) if max_date >= 100 else 12,
                    'end_day': int(str(max_date)[-2:]) if max_date >= 100 else 31
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

def create_filter_visualization(original_df: pd.DataFrame, input_df: pd.DataFrame, filtered_df: pd.DataFrame,
                                column: str, column_type: str):
    fig = go.Figure()

    base_count = len(filtered_df)
    input_count = len(input_df)

    if column_type == 'date':
        all_values = pd.to_numeric(input_df[column], errors='coerce').dropna()
        filtered_values = pd.to_numeric(filtered_df[column], errors='coerce').dropna()
        if all_values.empty:
            return html.Div("No valid date values found in the data")

        months = []
        for x in all_values:
            try:
                m = int(str(int(x))[:-2])  # mmdd → mm
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
                (float(f"{month}21"), float(f"{month}31")),
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
            barmode='overlay', showlegend=True
        )

    elif column_type == 'numeric':
        all_values = pd.to_numeric(input_df[column], errors='coerce')
        filtered_values = pd.to_numeric(filtered_df[column], errors='coerce')
        fig.add_trace(go.Histogram(x=all_values, name="data", opacity=0.5))
        fig.add_trace(go.Histogram(x=filtered_values, name="filterd data", opacity=0.7))
        fig.update_layout(barmode='overlay', showlegend=True)

    else:
        vc_all = input_df[column].value_counts()
        vc_fil = filtered_df[column].value_counts()
        fig.add_trace(go.Bar(x=vc_all.index.astype(str), y=vc_all.values, name='data', opacity=0.5))
        if not filtered_df.empty:
            fig.add_trace(go.Bar(x=vc_fil.index.astype(str), y=vc_fil.values, name='filterd data', opacity=0.7))
        fig.update_layout(showlegend=True)

    fig.update_layout(title=f"{column}", height=300, margin={'t': 30, 'b': 30, 'l': 30, 'r': 30})

    return html.Div([
        html.Strong(
            f"Filtered results: {base_count:,} out of total {input_count:,} cases",
            className="filter-result-text",
            style={
                'fontSize': '14px','color': '#2c3e50','display': 'block','marginBottom': '10px',
                'marginTop': '15px','padding': '8px','backgroundColor': '#f8f9fa',
                'borderRadius': '4px','border': '1px solid #dee2e6'
            }),
        dcc.Graph(figure=fig),
    ])

def apply_stored_filter(df: pd.DataFrame, filter_info: Dict[str, Any]) -> pd.DataFrame:
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

def create_filter_component(filter_id_num: int, column: str, column_type: str, df: pd.DataFrame):
    # 기존 UI 그대로 유지 (id: {'type': ..., 'index': filter_id_num})
    preview_data = df[column].head(10).tolist()  # 필요 시 사용
    months = [{'label': f'{i}', 'value': i} for i in range(1, 13)]
    days   = [{'label': f'{i}', 'value': i} for i in range(1, 32)]
    init   = get_initial_filter_values(df, column_type, column)

    # numeric 파라미터
    if column_type == 'numeric':
        rng = init.get('range', [0, 100])
        rmin, rmax = float(rng[0]), float(rng[1])
        rstep = (rmax - rmin) / 100 if rmin != rmax else 1
    else:
        rmin = 0; rmax = 100; rstep = 1

    label = f"Filter : {column}"
    return html.Div([
        html.Div([
            html.H4(label, style={'display': 'inline-block','fontSize': '16px','color': '#2c3e50','fontWeight': 'bold','margin': '0'}),
            html.Button('×', id={'type': 'remove-filter', 'index': filter_id_num},
                        style={'float': 'right','border': 'none','background': 'none',
                               'fontSize': '20px','cursor': 'pointer','color': '#dc3545','fontWeight': 'bold'})
        ], style={'marginBottom': '15px'}),

        # DATE
        html.Div([
            html.Div([
                html.Label('Start date:', style={'fontWeight': 'bold','marginRight': '20px','width': '80px','display': 'inline-block','verticalAlign': 'top','paddingTop': '5px','fontSize': '14px','color': '#495057'}),
                html.Div([
                    html.Div([
                        html.Label('Month:'), dcc.Dropdown(id={'type': 'filter-startdate-month','index': filter_id_num}, options=months, value=init.get('start_month'), placeholder='month', style={'width':'100px'})
                    ], style={'display': 'inline-block','marginRight': '10px'}),
                    html.Div([
                        html.Label('Day:'), dcc.Dropdown(id={'type': 'filter-startdate-day','index': filter_id_num}, options=days, value=init.get('start_day'), placeholder='day', style={'width':'100px'})
                    ], style={'display': 'inline-block'})
                ])
            ], style={'marginBottom': '10px'}),

            html.Div([
                html.Label('End date:', style={'fontWeight': 'bold','marginRight': '20px','width': '80px','display': 'inline-block','verticalAlign': 'top','paddingTop': '5px','fontSize': '14px','color': '#495057'}),
                html.Div([
                    html.Div([
                        html.Label('Month:'), dcc.Dropdown(id={'type': 'filter-enddate-month','index': filter_id_num}, options=months, value=init.get('end_month'), placeholder='month', style={'width':'100px'})
                    ], style={'display': 'inline-block','marginRight': '10px'}),
                    html.Div([
                        html.Label('Day:'), dcc.Dropdown(id={'type': 'filter-enddate-day','index': filter_id_num}, options=days, value=init.get('end_day'), placeholder='day', style={'width':'100px'})
                    ], style={'display': 'inline-block'})
                ])
            ])
        ], id={'type': 'date-container', 'index': filter_id_num}, style={'display': 'block' if column_type == 'date' else 'none'}),

        # NUMERIC
        html.Div([
            dcc.RangeSlider(
                id={'type': 'numeric-range', 'index': filter_id_num},
                min=rmin, max=rmax, step=rstep, value=[rmin, rmax],
                tooltip={'placement': 'bottom', 'always_visible': True},
                marks={rmin: f'{rmin:.2f}', rmax: f'{rmax:.2f}'}
            )
        ], id={'type': 'numeric-container','index': filter_id_num}, style={'display': 'block' if column_type == 'numeric' else 'none'}),

        # CATEGORICAL
        html.Div([
            dcc.Dropdown(
                id={'type': 'categorical-filter', 'index': filter_id_num},
                options=[{'label': str(v), 'value': str(v)} for v in sorted(df[column].dropna().unique())],
                multi=True,
                placeholder="Select categories...",
                value=[str(v) for v in sorted(df[column].dropna().unique())]
            )
        ], id={'type': 'categorical-container','index': filter_id_num}, style={'display': 'block' if column_type == 'categorical' else 'none'}),

        html.Div(id={'type': 'filter-result','index': filter_id_num}),
        html.Hr()
    ], className="filter-block")
