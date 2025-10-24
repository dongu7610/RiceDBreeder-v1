"""
GWAS Module with Enhanced Scatter Plot Features

This module provides GWAS (Genome-Wide Association Studies) analysis functionality
with enhanced scatter plot visualization capabilities integrated from gwas_module2.py.

Key Features:
1. Original sample-based scatter plots (Y-axis: sample names)
2. New P-value scatter plots (Y-axis: -log10(p-value))
3. Unique position filtering for cross-sample comparison
4. P-value cutoff visualization
5. MAF (Minor Allele Frequency) filtering
6. Interactive trait selection and highlighting

Usage Example:
    # Create GWAS module instance
    gwas_module = GWASModule()
    
    # Bar chart (original functionality)
    fig_bar = gwas_module.make_figure_for_subtrait(subtrait, selected_samples, hide_vals, selected_traits)
    
    # P-value scatter plot (NEW - from gwas_module2.py)
    fig_scatter = gwas_module.make_pvalue_scatter_for_samples(
        selected_samples=['sample1', 'sample2'],
        selected_subtraits=['trait1', 'trait2'],
        maf_threshold=0.05,
        input_loci=['chr1:12345', 'chr2:67890'],
        pvalue_cutoff=0.01,
        unique_only=False  # Set True for sample-unique positions only
    )
    
    # Original sample scatter plot (Y-axis: sample names)
    fig_original = gwas_module.make_scatter_for_samples(
        selected_samples=['sample1', 'sample2'],
        selected_subtraits=['trait1', 'trait2'],
        maf_threshold=0.05,
        input_loci=['chr1:12345', 'chr2:67890']
    )

Column Requirements:
    Sample DataFrames must contain:
    - Chromosome: Chromosome identifier (str)
    - Start_Position/Position: Genomic position (int)
    - Subtrait: Trait category (str)
    - P-value: Statistical significance (float)
    - MAF: Minor Allele Frequency (float)
    - Trait: Specific trait name (str)
    - Locus: Combined "chr:pos" format (str)
"""

import os
import glob
import json
import pandas as pd
import dash
from dash import dcc, html, dash_table, callback_context, no_update
from dash.dependencies import Input, Output, State, ALL
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import plotly.subplots as sp
import plotly.express as px
import numpy as np


# 마커 심볼 정의 (gwas_module2.py에서 가져옴)
MARKER_SYMBOLS = [
    'circle', 'square', 'diamond', 'cross', 'x', 'triangle-up', 'triangle-down',
    'triangle-left', 'triangle-right', 'star', 'pentagon', 'hexagon', 'hexagon2', 'octagon', 'plus'
]

ALL_SUBTRAITS = [
    'Yield(14)',  'Stress(66)', 'Plantvigor(23)', 'Biochemical(103)', 'Plantgrowth_Development(17)', 'Plant_quality(26)', 'Biologicalprocess(19)', 'Plantmorphology(191)', 'Sterility_Fertility(2)'
]
ALL_SUBTRAIT_COLORS = px.colors.qualitative.Set1[:9]
SUBTRAIT_COLOR_MAP = {sub: color for sub, color in zip(ALL_SUBTRAITS, ALL_SUBTRAIT_COLORS)}

class GWASModule:
    def __init__(self, data_path=None ,filtered_results_path='./data_genotype/'):
        """
        Initialize GWAS Module
        
        Args:
            data_path: Path to the main GWAS CSV file
            filtered_results_path: Path to filtered results directory
        """
        self.data_path = data_path
        self.filtered_results_path = filtered_results_path
        self.full_df = None
        self.sample_dfs = {}
        self.samples = []
        self.color_map = {}
        self.app = None
        
        # Load data during initialization
        
        self._load_data()
        self._setup_colors()
        #self._index_samples()
    '''def _index_samples(self):
        tsv_paths = glob.glob(os.path.join(self.filtered_results_path, '*_SNP_gwas_yes.tsv'))
        names = []
        for path in tsv_paths:
            name = os.path.basename(path).replace('_SNP_gwas_yes.tsv', '')
            names.append(name)
        self.samples = sorted(set(names))

    def get_samples(self):
        return list(self.samples)

    def get_sample_data(self, sample_name: str) -> pd.DataFrame:
        """요청 시점에 해당 샘플 TSV만 읽어서 반환(지연 로딩)"""
        if sample_name in self.sample_dfs:
            return self.sample_dfs[sample_name]

        path = os.path.join(self.filtered_results_path, f'{sample_name}_SNP_gwas_yes.tsv')
        if not os.path.exists(path):
            return pd.DataFrame(columns=[
                'Variant ID','Chromosome','Location','Allele','MAF','Minor Allele',
                'Trait','P-value','R2 (%)','Lead SNP','PMID','Genes','Subtrait',
                'Parsed_Genes','Start_Position','sample_check1','sample_check2',
                'Locus','Position'
            ])

        tmp = pd.read_csv(
            path, sep='\t', header=None,
            names=[
                'index','Variant ID','Chromosome','Location','Allele','MAF',
                'Minor Allele','Trait','P-value','R2 (%)','Lead SNP','PMID',
                'Genes','Subtrait','Parsed_Genes','Start_Position',
                'sample_check1','sample_check2'
            ],
            dtype=str
        )
        tmp['Locus'] = tmp['Chromosome'].astype(str) + ':' + tmp['Start_Position'].astype(str)
        tmp['Position'] = tmp['Start_Position'].astype(str).astype(int, errors='ignore')
        # 필요 시 전처리(숫자형 변환)
        for col in ['MAF', 'P-value']:
            if col in tmp.columns:
                tmp[col] = pd.to_numeric(tmp[col], errors='coerce')
        self.sample_dfs[sample_name] = tmp
        return tmp
        '''
    
    def _load_data(self):
        """Load GWAS data and sample data"""
        # Load main GWAS data
        '''usecols = ['Chromosome', 'Start_Position', 'Trait', 'Subtrait', 'Minor Allele', 'MAF']
        try:
            self.full_df = pd.read_csv(
                self.data_path,
                usecols=usecols + ['Variant ID', 'P-value', 'PMID'],
                dtype=str
            )
            self.full_df['Locus'] = self.full_df['Chromosome'] + ':' + self.full_df['Start_Position']
            self.full_df['Position'] = self.full_df['Start_Position'].astype(int)
        except FileNotFoundError:
            self.full_df = pd.DataFrame(columns=usecols + ['Variant ID', 'P-value', 'PMID', 'Locus', 'Position'])
        '''
        # Load per-sample GWAS filtered data
        tsv_paths = glob.glob(f'{self.filtered_results_path}*_SNP_gwas_yes.tsv')
        for path in tsv_paths:
            name = os.path.basename(path).replace('_SNP_gwas_yes.tsv', '')
            tmp = pd.read_csv(
                path, sep='\t', header=None,
                names=[
                    'index', 'Variant ID', 'Chromosome', 'Location', 'Allele', 'MAF',
                    'Minor Allele', 'Trait', 'P-value', 'R2 (%)', 'Lead SNP', 'PMID',
                    'Genes', 'Subtrait', 'Parsed_Genes', 'Start_Position',
                    'sample_check1', 'sample_check2'
                ],
                dtype=str
            )
            tmp['Locus'] = tmp['Chromosome'] + ':' + tmp['Start_Position']
            tmp['Position'] = tmp['Start_Position'].astype(int)
            self.sample_dfs[name] = tmp

        self.samples = sorted(self.sample_dfs.keys())
    

    def _setup_colors(self):
        """Setup color mapping for samples and GWAS"""
        color_palette = px.colors.qualitative.Plotly
        self.color_map = {"GWAS": "black"}
        for i, sample in enumerate(self.samples):
            self.color_map[sample] = color_palette[i % len(color_palette)]
    
    def create_app(self, external_stylesheets=None):
        """
        Create and configure the Dash application
        
        Args:
            external_stylesheets: List of external stylesheets
            
        Returns:
            dash.Dash: Configured Dash application
        """
        if external_stylesheets is None:
            external_stylesheets = [dbc.themes.BOOTSTRAP]
            
        self.app = dash.Dash(__name__, external_stylesheets=external_stylesheets, suppress_callback_exceptions=True)
        self.app.layout = self._create_layout()
        self._register_callbacks()
        
        return self.app
    
    def _create_layout(self):
        """Create the dashboard layout"""
        return html.Div([
            html.H2("GWAS Sample Comparison Dashboard"),

            # 1. Sample selector
            html.Div([
                html.Label("Select Samples:"),
                dcc.Dropdown(
                    id='sample-dropdown',
                    options=[{'label': s, 'value': s} for s in self.samples],
                    multi=True,
                    placeholder='Choose one or more samples'
                )
            ], style={'width': '40%', 'margin-bottom': '20px'}),
            
            # 2. Subtrait dropdown (now multi-select)
            html.Div([
                html.Label("Choose Subtraits:"),
                dcc.Dropdown(
                    id='subtrait-dropdown', 
                    placeholder='Select one or more subtraits',
                    multi=True
                )
            ], style={'width': '40%', 'margin-bottom': '20px'}),
            
            html.Hr(),
            
            # 3. Trait bar chart (with toggle visibility and icon) - initially hidden
            html.Div([
                dbc.Button([
                    html.Span("▼", id="trait-chart-icon", style={'margin-right': '10px'}),
                    "Trait Counts Chart"
                ], 
                    id="toggle-trait-chart",
                    color="link",
                    style={'padding': '10px 0', 'text-align': 'left', 'border': 'none'}
                ),
                dbc.Collapse(
                    id="trait-chart-collapse",
                    is_open=True,
                    children=[
                        # Hide traits checkbox - only visible when trait chart is open
                        html.Div([
                            dcc.Checklist(
                                id='hide-traits-checkbox',
                                options=[{'label': 'Hide traits not in selected samples', 'value': 'hide'}],
                                value=[],
                                style={'margin-left': '20px', 'margin-bottom': '10px'}
                            )
                        ], id='hide-checkbox-container'),
                        html.Div(
                            id='trait-bar-container',
                            children=[],
                            style={'display': 'none'}  # Initially hidden
                        )
                    ]
                )
            ]),
            
            # 4. Selected traits store and display
            dcc.Store(id='selected-traits-store', data=[]),
            html.Div(id='selected-traits-container', style={'margin': '10px 0'}),
            
            # 5. Sample scatter (with toggle visibility and icon) - initially hidden
            html.Div([
                dbc.Button([
                    html.Span("▼", id="scatter-icon", style={'margin-right': '10px'}),
                    "Sample Variant Positions"
                ], 
                    id="toggle-scatter",
                    color="link",
                    style={'padding': '10px 0', 'text-align': 'left', 'border': 'none'}
                ),
                dbc.Collapse(
                    id="scatter-collapse",
                    is_open=True,
                    children=[
                        html.Div(
                            id='sample-scatter-container',
                            children=[dcc.Graph(id='sample-scatter')],
                            style={'display': 'none'}  # Initially hidden
                        )
                    ]
                )
            ]),
            
            html.Hr(),
            
            # 6. Data table and download section
            html.Div([
                html.H4("Selected Loci Data"),
                # Download button - only visible when samples are selected
                html.Div([
                    dbc.Button("Download All Loci with Samples", id="download-all-btn", color="primary", style={'margin-bottom': '10px', 'margin-right': '10px'}),
                    dbc.Button("Copy All Loci", id="copy-loci-btn", color="secondary", style={'margin-bottom': '10px'}),
                    dcc.Download(id="download-all-loci")
                ], id='download-container', style={'display': 'none'}),
                
                # Copy result display
                html.Div(id='copy-result', style={'margin-bottom': '10px'}),
                
                # Store component for filter state management
                dcc.Store(id='filter-states', data={'maf': None}),
                dcc.Store(id='input-loci-store', data=[]),
                
                # Filter section
                html.Div([
                    # MAF filter button (initially displayed)
                    html.Button("Consider MAF", id="maf-consider-btn", n_clicks=0, 
                               style={'margin-right': '10px'}),
                    
                    # MAF input field (initially hidden, dynamically added/removed)
                    html.Div(id="maf-input-container", children=[], style={'display': 'inline-block'}),
                    
                    # MAF span (initially hidden, dynamically added/removed)
                    html.Div(id="maf-span-container", children=[], style={'display': 'inline-block'}),
                    
                    # Locus input button (initially displayed)
                    html.Button("Input Locus", id="locus-input-btn", n_clicks=0, 
                               style={'margin-right': '10px', 'margin-left': '20px'}),
                    
                    # Locus input field (initially hidden, dynamically added/removed)
                    html.Div(id="locus-input-container", children=[], style={'display': 'inline-block'}),
                    
                    # Locus span (initially hidden, dynamically added/removed)
                    html.Div(id="locus-span-container", children=[], style={'display': 'inline-block'}),
                    
                ], style={'margin-top': '10px', 'margin-bottom': '10px'}),
                
                dash_table.DataTable(
                    id='locus-table',
                    columns=[],
                    data=[],
                    page_size=20,
                    style_table={'overflowX': 'auto'},
                    style_cell={'textAlign': 'left', 'padding': '10px'},
                    style_header={'backgroundColor': 'lightgrey', 'fontWeight': 'bold'},
                    style_data_conditional=[
                        {
                            'if': {'column_id': 'inputlocus', 'filter_query': '{inputlocus} = Yes'},
                            'backgroundColor': '#e8f4fd',
                            'color': 'black',
                        },
                        {
                            'if': {'column_id': 'Source', 'filter_query': '{Source} = inputlocus'},
                            'backgroundColor': '#e8f4fd',
                            'color': 'black',
                        }
                    ]
                )
            ], style={'margin-top': '20px'})
        ])

    def run_server(self, debug=True, port=8050, host='0.0.0.0'):
        """
        Run the Dash server
        
        Args:
            debug: Enable debug mode
            port: Port number
            host: Host address
        """
        if self.app is None:
            self.create_app()
        
        self.app.run_server(debug=debug, port=port, host=host)
    
    def get_samples(self):
        """Get list of available samples"""
        return self.samples.copy()
    
    def get_sample_data(self, sample_name):
        """Get data for a specific sample"""
        return self.sample_dfs.get(sample_name, pd.DataFrame())
    
    def get_gwas_data(self):
        """Get main GWAS data"""
        return self.full_df.copy()
    
    def filter_unique_positions(self, data, selected_samples):
        """
        Filter data to show only unique positions per sample
        
        Args:
            data: DataFrame with GWAS data
            selected_samples: List of selected sample names
            
        Returns:
            DataFrame: Filtered data with unique positions
        """
        if 'Sample' not in data.columns:
            return data
            
        if len(selected_samples) > 1:
            # 여러 샘플: 각 position별로 샘플 수 계산하여 unique position 찾기
            position_cols = ['Chromosome', 'Start_Position'] if 'Start_Position' in data.columns else ['Chromosome', 'Position']
            if all(col in data.columns for col in position_cols):
                position_counts = data.groupby(position_cols)['Sample'].nunique().reset_index()
                position_counts.columns = position_cols + ['sample_count']
                
                # 하나의 샘플에만 있는 position 찾기
                unique_positions = position_counts[position_counts['sample_count'] == 1][position_cols]
                
                # unique position만 필터링
                filtered_data = data.merge(unique_positions, on=position_cols, how='inner')
                return filtered_data
        else:
            # 단일 샘플: 중복 position 제거 (가장 유의한 P-value만 유지)
            position_cols = ['Chromosome', 'Start_Position'] if 'Start_Position' in data.columns else ['Chromosome', 'Position']
            if 'P-value' in data.columns and all(col in data.columns for col in position_cols):
                # P-value를 숫자로 변환
                data['P-value'] = pd.to_numeric(data['P-value'], errors='coerce')
                # 각 position에서 가장 작은 P-value만 유지
                filtered_data = data.loc[data.groupby(position_cols)['P-value'].idxmin()]
                return filtered_data
        
        return data
    
    def _register_callbacks(self):
        """Register all Dash callbacks"""
        
        # Toggle trait chart visibility with icon
        @self.app.callback(
            [Output("trait-chart-collapse", "is_open"),
             Output("trait-chart-icon", "children")],
            Input("toggle-trait-chart", "n_clicks"),
            State("trait-chart-collapse", "is_open")
        )
        def toggle_trait_chart(n_clicks, is_open):
            if n_clicks:
                new_state = not is_open
                icon = "▼" if new_state else "▶"
                return new_state, icon
            return is_open, "▼" if is_open else "▶"

        # Control chart visibility based on selections
        @self.app.callback(
            Output('trait-bar-container', 'style'),
            [Input('sample-dropdown', 'value'),
             Input('subtrait-dropdown', 'value')]
        )
        def control_trait_chart_visibility(selected_samples, selected_subtraits):
            # Show trait chart only when both samples and subtraits are selected
            if selected_samples and selected_subtraits:
                return {'display': 'block'}
            else:
                return {'display': 'none'}

        # Unified scatter plot callback
        @self.app.callback(
            [
                Output('sample-scatter', 'figure'),
                Output('sample-scatter-container', 'style'),
                Output('scatter-collapse', 'is_open'),
                Output("scatter-icon", "children")
            ],
            [
                Input('sample-dropdown', 'value'),
                Input('subtrait-dropdown', 'value'),
                Input('selected-traits-store', 'data'),
                Input("toggle-scatter", "n_clicks"),
                Input('filter-states', 'data'),
                Input('input-loci-store', 'data')
            ],
            [State("scatter-collapse", "is_open")]
        )
        def update_scatter_unified(selected_samples, selected_subtraits, selected_traits, toggle_clicks, filter_states, input_loci, current_is_open):
            ctx = callback_context
            
            # MAF 필터 값 추출
            maf_threshold = filter_states.get('maf')
            
            # Check if toggle button was clicked
            if ctx.triggered and ctx.triggered[0]['prop_id'] == 'toggle-scatter.n_clicks':
                # Manual toggle - just flip the current state
                new_is_open = not current_is_open
                icon = "▼" if new_is_open else "▶"
                
                # Keep current figure and style based on selections
                if not selected_samples and not input_loci:
                    return go.Figure(), {'display': 'none'}, new_is_open, icon
                else:
                    fig = self.make_scatter_for_samples(selected_samples, selected_subtraits, maf_threshold, input_loci)
                    self.highlight_traits_on_scatter(fig, selected_traits)
                    return fig, {'display': 'block'}, new_is_open, icon
            
            # Automatic control based on selections
            # 1) 샘플과 input loci 모두 없으면 숨김
            if not selected_samples and not input_loci:
                return go.Figure(), {'display': 'none'}, False, "▶"

            # 2) 샘플이나 input loci가 있으면 표시
            fig = self.make_scatter_for_samples(selected_samples, selected_subtraits, maf_threshold, input_loci)
            self.highlight_traits_on_scatter(fig, selected_traits)
            return fig, {'display': 'block'}, True, "▼"

        # Populate subtrait dropdown
        @self.app.callback(
            Output('subtrait-dropdown', 'options'),
            Input('sample-dropdown', 'value')
        )
        def update_subtrait_options(_):
            st_list = self.full_df['Subtrait'].dropna().astype(str).unique().tolist()
            st_list.sort()
            return [{'label': st, 'value': st} for st in st_list]

        # Generate individual trait-bar-charts for each subtrait
        @self.app.callback(
            Output('trait-bar-container', 'children'),
            [Input('sample-dropdown', 'value'),
             Input('subtrait-dropdown', 'value'),
             Input('hide-traits-checkbox', 'value'),
             Input('selected-traits-store', 'data')]
        )
        def update_trait_bar_container(selected_samples, selected_subtraits, hide_vals, selected_traits):
            if not selected_subtraits:
                return []
            
            charts = []
            for subtrait in selected_subtraits:
                height = self.calc_height(subtrait, selected_samples, hide_vals, selected_traits)
                figure = self.make_figure_for_subtrait(subtrait, selected_samples, hide_vals, selected_traits)
                
                charts.append(
                    dcc.Graph(
                        id={'type': 'trait-bar', 'subtrait': subtrait},
                        figure=figure,
                        style={'height': f'{height}px', 'margin-bottom': '20px'}
                    )
                )
            
            return charts

        # Handle trait selection from individual charts
        @self.app.callback(
            Output('selected-traits-store', 'data'),
            [Input({'type': 'trait-bar', 'subtrait': ALL}, 'clickData'),
             Input({'type': 'remove-button', 'trait': ALL}, 'n_clicks')],
            [State('selected-traits-store', 'data')]
        )
        def handle_trait_selection(click_data_list, remove_clicks, current_selected):
            ctx = callback_context
            if not ctx.triggered:
                return current_selected or []
            
            trigger = ctx.triggered[0]['prop_id']
            
            # Handle remove button clicks
            if 'remove-button' in trigger:
                trait = json.loads(trigger.split('.')[0])['trait']
                return [t for t in (current_selected or []) if t != trait]
            
            # Handle chart clicks
            if click_data_list and any(click_data_list):
                for click_data in click_data_list:
                    if click_data and 'points' in click_data:
                        for point in click_data['points']:
                            if 'y' in point:
                                clicked_trait = point['y']
                                current_selected = current_selected or []
                                
                                if clicked_trait in current_selected:
                                    current_selected.remove(clicked_trait)
                                else:
                                    current_selected.append(clicked_trait)
                                
                                return current_selected
            
            return current_selected or []

        # Render selected traits
        @self.app.callback(
            Output('selected-traits-container', 'children'),
            Input('selected-traits-store', 'data')
        )
        def render_selected_traits(selected_traits):
            spans = []
            for trait in selected_traits:
                spans.append(html.Span([
                    trait,
                    html.Button('×', id={'type': 'remove-button', 'trait': trait}, n_clicks=0,
                                style={'margin-left': '5px', 'background': 'transparent', 'border': 'none'})
                ], style={
                    'display': 'inline-block', 'padding': '5px', 'border': '1px solid #ccc',
                    'border-radius': '4px', 'margin-right': '5px', 'margin-bottom': '5px'
                }))
            return spans

        # Data table callback
        @self.app.callback(
            [Output('locus-table', 'columns'),
             Output('locus-table', 'data')],
            [Input('selected-traits-store', 'data'),
             Input('sample-dropdown', 'value'),
             Input('subtrait-dropdown', 'value'),
             Input('filter-states', 'data'),
             Input('input-loci-store', 'data')]
        )
        def update_locus_table(selected_traits, selected_samples, selected_subtraits, filter_states, input_loci):
            return self.update_locus_table_data(selected_traits, selected_samples, selected_subtraits, filter_states, input_loci)

        # Add other callbacks (download, copy, MAF filter, locus input, etc.)
        self._register_additional_callbacks()
    
    def make_scatter_for_samples(self, selected_samples, selected_subtraits=None, maf_threshold=None, input_loci=None):
        """Generate scatter plot for selected samples and input loci"""
        # Create 4×3 grid for chromosomes
        fig = sp.make_subplots(
            rows=3, cols=4,
            subplot_titles=[f'chr{i}' for i in range(1, 13)],
            shared_xaxes='columns',
            shared_yaxes='rows',
            vertical_spacing=0.1,
            horizontal_spacing=0.05
        )
        
        # Define y-axis categories (inputlocus + samples)
        y_categories = []
        if input_loci:
            y_categories.append('inputlocus')
        if selected_samples:
            y_categories.extend(selected_samples)
        
        # Process each chromosome
        for chrom_idx in range(1, 13):
            row = ((chrom_idx-1) // 4) + 1
            col = ((chrom_idx-1) % 4) + 1
            chrom_str = str(chrom_idx)
            
            # Add input loci data first
            if input_loci:
                input_positions = []
                input_traits = []
                
                for locus in input_loci:
                    # chr 접두사 제거
                    locus_original = locus.replace('chr', '') if locus.startswith('chr') else locus
                    
                    # chromosome과 position 분리
                    if ':' in locus_original:
                        locus_chrom, locus_pos = locus_original.split(':')
                        if locus_chrom == chrom_str:
                            # 샘플에서 먼저 찾기
                            trait_found = None
                            if selected_samples:
                                for sample in selected_samples:
                                    sample_data = self.sample_dfs[sample][self.sample_dfs[sample]['Locus'] == locus_original]
                                    if not sample_data.empty:
                                        trait_found = sample_data.iloc[0]['Trait']
                                        break
                            
                            # 샘플에 없으면 GWAS에서 찾기
                            if not trait_found:
                                gwas_data = self.full_df[self.full_df['Locus'] == locus_original]
                                if not gwas_data.empty:
                                    trait_found = gwas_data.iloc[0]['Trait']
                            
                            input_positions.append(int(locus_pos))
                            input_traits.append(trait_found or 'Unknown')
                
                if input_positions:
                    fig.add_trace(
                        go.Scattergl(
                            x=input_positions,
                            y=['inputlocus'] * len(input_positions),
                            mode='markers',
                            name='inputlocus',
                            marker=dict(
                                color='red',
                                size=10, 
                                line=dict(width=2),
                                opacity=1.0,
                                symbol='diamond'
                            ),
                            customdata=input_traits,
                            showlegend=(chrom_idx == 1)
                        ),
                        row=row, col=col
                    )
            
            # Add sample data
            if selected_samples:
                for sample in selected_samples:
                    df_s = self.sample_dfs[sample]
                    # Get loci for this chromosome
                    sample_df_filtered = df_s[df_s['Chromosome'] == chrom_str]
                    
                    # Apply MAF filter if exists
                    if maf_threshold is not None:
                        sample_df_filtered = sample_df_filtered[sample_df_filtered['MAF'].astype(float) > maf_threshold]
                    
                    positions = sorted(sample_df_filtered['Position'].astype(int).tolist())
                    sample_traits = sample_df_filtered['Trait'].tolist()
                    
                    if positions:
                        fig.add_trace(
                            go.Scattergl(
                                x=positions,
                                y=[sample] * len(positions),
                                mode='markers',
                                name=sample,
                                marker=dict(
                                    color=self.color_map[sample], 
                                    size=8, 
                                    line=dict(width=0.1),
                                    opacity=1.0
                                ),
                                customdata=sample_traits,
                                showlegend=(chrom_idx == 1)
                            ),
                            row=row, col=col
                        )
        
        # Update layout
        fig.update_layout(
            height=600,
            showlegend=True,
            title="Sample Variant Positions by Chromosome – Highlight Selected Traits"
        )
        
        # Set y-axis categories for all subplots
        for i in range(1, 13):
            row = ((i-1) // 4) + 1
            col = ((i-1) % 4) + 1
            fig.update_yaxes(categoryorder='array', categoryarray=y_categories, row=row, col=col)
        
        return fig

    def highlight_traits_on_scatter(self, fig, selected_traits):
        """Highlight selected traits on the scatter plot by changing line width only"""
        # Reset all to default line width first
        for trace in fig.data:
            if hasattr(trace, 'marker') and trace.marker:
                trace.marker.line.width = 1
                trace.marker.opacity = 1.0
        
        if not selected_traits:
            return
        
        # Highlight selected traits with thicker line width
        for trace in fig.data:
            if hasattr(trace, 'marker') and trace.marker and hasattr(trace, 'customdata'):
                if trace.customdata:
                    # Create line width array for each point
                    line_widths = []
                    for i, trait in enumerate(trace.customdata):
                        if trait in selected_traits:
                            line_widths.append(1)
                        else:
                            line_widths.append(0)
                    
                    trace.marker.line.width = line_widths 

    def calc_height(self, subtrait, selected_samples, hide_vals, selected_traits):
        """Calculate height for individual subtrait chart"""
        MIN_HEIGHT = 300
        MAX_HEIGHT = 800
        BAR_HEIGHT = 20
        PADDING = 100
        
        # Get all traits for this subtrait
        df_sub = self.full_df[self.full_df['Subtrait'] == subtrait]
        all_traits = sorted(df_sub['Trait'].unique().tolist())
        
        if 'hide' in hide_vals and selected_samples:
            keep = []
            for trait in all_traits:
                present = False
                for sample in selected_samples:
                    if not self.sample_dfs[sample][
                        (self.sample_dfs[sample]['Subtrait'] == subtrait) & 
                        (self.sample_dfs[sample]['Trait'] == trait)
                    ].empty:
                        present = True
                        break
                if present:
                    keep.append(trait)
            traits_list = keep
        else:
            traits_list = all_traits
        
        n = len(traits_list)
        h = max(MIN_HEIGHT, min(n * BAR_HEIGHT + PADDING, MAX_HEIGHT))
        return h

    def make_figure_for_subtrait(self, subtrait, selected_samples, hide_vals, selected_traits):
        """Create individual figure for one subtrait"""
        BAR_THICKNESS = 0.6
        
        fig = go.Figure()
        
        # Get all traits for this subtrait
        df_sub = self.full_df[self.full_df['Subtrait'] == subtrait]
        all_traits = sorted(df_sub['Trait'].unique().tolist())
        
        if 'hide' in hide_vals and selected_samples:
            keep = []
            for trait in all_traits:
                present = False
                for sample in selected_samples:
                    if not self.sample_dfs[sample][
                        (self.sample_dfs[sample]['Subtrait'] == subtrait) & 
                        (self.sample_dfs[sample]['Trait'] == trait)
                    ].empty:
                        present = True
                        break
                if present:
                    keep.append(trait)
            traits_list = keep
        else:
            traits_list = all_traits
        
        if not traits_list:
            return fig
        
        # Calculate trait union counts for samples
        trait_union_counts = {}
        for trait in traits_list:
            union_loci = set()
            for sample in selected_samples or []:
                sample_loci = set(self.sample_dfs[sample][
                    (self.sample_dfs[sample]['Subtrait'] == subtrait) & 
                    (self.sample_dfs[sample]['Trait'] == trait)
                ]['Locus'].tolist())
                union_loci.update(sample_loci)
            trait_union_counts[trait] = len(union_loci)
        
        # GWAS counts
        gwas_counts = [len(df_sub[df_sub['Trait'] == trait]['Locus'].unique()) for trait in traits_list]
        gwas_colors = ['darkred' if trait in selected_traits else 'lightcoral' for trait in traits_list]
        
        # Add GWAS bars with controlled thickness
        fig.add_trace(go.Bar(
            y=traits_list,
            x=gwas_counts,
            orientation='h',
            name="GWAS",
            marker=dict(
                color=gwas_colors,
                line=dict(width=1)
            ),
            width=BAR_THICKNESS,
            showlegend=True
        ))
        
        # Add sample union contribution bars (stacked) with highlighting
        if selected_samples:
            union_counts = [trait_union_counts[trait] for trait in traits_list]
            union_colors = ['darkred' if trait in selected_traits else 'steelblue' for trait in traits_list]
            
            fig.add_trace(go.Bar(
                y=traits_list,
                x=union_counts,
                orientation='h',
                name="Sample Union",
                marker=dict(
                    color=union_colors,
                    line=dict(width=1)
                ),
                width=BAR_THICKNESS,
                showlegend=True
            ))
        
        fig.update_layout(
            barmode='stack',
            clickmode='event+select',
            title=f"Trait Counts for Subtrait: {subtrait}",
            margin=dict(l=200, r=50, t=80, b=50),
            bargap=0.2,
            bargroupgap=0.1
        )
        
        return fig

    def update_locus_table_data(self, selected_traits, selected_samples, selected_subtraits, filter_states, input_loci):
        """Update locus table data"""
        if not selected_samples:
            return [], []
        
        # MAF 필터 값 추출
        maf_threshold = filter_states.get('maf')
        
        # Collect all sample loci
        all_loci = set()
        sample_data_map = {}
        
        # 수정: 선택된 traits가 있든 없든 샘플의 모든 데이터 포함
        for sample in selected_samples:
            df_s = self.sample_dfs[sample]
            
            # trait 필터링: 선택된 traits가 있으면 필터링, 없으면 전체 포함
            if selected_traits:
                df_filt = df_s[df_s['Trait'].isin(selected_traits)]
            else:
                df_filt = df_s.copy()  # 전체 샘플 데이터 포함
            
            # Apply MAF filter if exists (교집합)
            if maf_threshold is not None:
                df_filt = df_filt[df_filt['MAF'].astype(float) > maf_threshold]
            
            # Remove duplicates by Locus to ensure unique index
            df_filt = df_filt.drop_duplicates(subset='Locus', keep='first')
            
            sample_loci = set(df_filt['Locus'].tolist())
            all_loci.update(sample_loci)
            
            # Store sample data for lookup - now with unique Locus index
            sample_data_map[sample] = df_filt.set_index('Locus').to_dict('index')
        
        # Input loci 추가 처리
        input_loci_processed = []
        if input_loci:
            for locus in input_loci:
                # chr 접두사 제거해서 원본 형태로 변환 (테이블에서는 chr 접두사 붙여서 표시)
                locus_original = locus.replace('chr', '') if locus.startswith('chr') else locus
                input_loci_processed.append(locus_original)
                all_loci.add(locus_original)
        
        if not all_loci:
            return [], []
        
        # Build table rows for each unique locus
        rows = []
        for locus in all_loci:
            is_input_locus = locus in input_loci_processed
            
            # 기존 샘플 데이터에서 locus 정보 찾기
            locus_info = None
            for sample in selected_samples:
                if locus in sample_data_map[sample]:
                    locus_info = sample_data_map[sample][locus]
                    break
            
            # 샘플에 없으면 GWAS에서 가져오기 (input loci의 경우)
            if not locus_info and is_input_locus:
                gwas_data = self.full_df[self.full_df['Locus'] == locus]
                if not gwas_data.empty:
                    locus_info = gwas_data.iloc[0].to_dict()
            
            if not locus_info:
                continue
                
            row_data = {
                'Source': 'inputlocus' if is_input_locus else 'Sample',
                'inputlocus': 'Yes' if is_input_locus else 'No',
                'Variant ID': locus_info.get('Variant ID', ''),
                'Locus': f"chr{locus}",
                'Minor Allele': locus_info.get('Minor Allele', ''),
                'MAF': locus_info.get('MAF', ''),
                'P-value': locus_info.get('P-value', ''),
                'PMID': locus_info.get('PMID', ''),
                'Trait': locus_info.get('Trait', ''),
                'Subtrait': locus_info.get('Subtrait', '')
            }
            
            # Add sample presence and value columns
            for sample in selected_samples:
                if locus in sample_data_map[sample]:
                    #row_data[f'{sample}_present'] = sample
                    row_data[f'{sample}_value'] = sample_data_map[sample][locus].get('sample_check2', '')
                else:
                    # 샘플에 없는 경우 빈 값
                    #row_data[f'{sample}_present'] = ""
                    row_data[f'{sample}_value'] = ""
            
            rows.append(row_data)
        
        if not rows:
            return [], []
        
        # Convert to DataFrame for sorting
        table_df = pd.DataFrame(rows)
        
        # Chromosome sorting logic - extract chromosome number from Locus
        table_df['chr_order'] = (
            table_df['Locus']
              .str.replace(r'^chr', '', regex=True)
              .str.split(':').str[0]
              .astype(int)
        )
        # Extract position for sorting
        table_df['pos_order'] = (
            table_df['Locus']
              .str.split(':').str[1]
              .astype(int)
        )
        table_df = table_df.sort_values(['chr_order', 'pos_order'])
        table_df = table_df.drop(columns=['chr_order', 'pos_order'])
        
        # Convert back to records
        rows = table_df.to_dict('records')
        
        # Create columns in the requested order
        base_columns = ['Source', 'inputlocus', 'Variant ID', 'Locus', 'Minor Allele', 'MAF', 'P-value', 'PMID', 'Trait', 'Subtrait']
        sample_columns = []
        for sample in selected_samples:
            sample_columns.extend([f'{sample}_value'])#f'{sample}_present',
        
        all_columns = base_columns + sample_columns
        
        columns = [{'name': col, 'id': col} for col in all_columns]
        
        return columns, rows

    def _register_additional_callbacks(self):
        """Register additional callbacks for download, copy, MAF filter, etc."""
        
        # Download all sample loci callback
        @self.app.callback(
            Output("download-all-loci", "data"),
            Input("download-all-btn", "n_clicks"),
            [State('sample-dropdown', 'value'),
             State('subtrait-dropdown', 'value'),
             State('filter-states', 'data'),
             State('input-loci-store', 'data')],
            prevent_initial_call=True
        )
        def download_all_loci(n_clicks, selected_samples, selected_subtraits, filter_states, input_loci):
            if not selected_samples and not input_loci:
                return no_update
            
            # Implementation similar to original but using self methods
            return self._handle_download_all_loci(selected_samples, selected_subtraits, filter_states, input_loci)

        # Control hide-checkbox visibility based on trait chart collapse state


        # Control download button visibility based on selected samples or input loci
        @self.app.callback(
            Output('download-container', 'style'),
            [Input('sample-dropdown', 'value'),
             Input('input-loci-store', 'data')]
        )
        def control_download_visibility(selected_samples, input_loci):
            if (selected_samples and len(selected_samples) > 0) or (input_loci and len(input_loci) > 0):
                return {'margin-bottom': '10px', 'display': 'block'}
            else:
                return {'display': 'none'}

        # Copy selected loci callback
        @self.app.callback(
            Output('copy-result', 'children'),
            [Input('copy-loci-btn', 'n_clicks')],
            [State('locus-table', 'data')],
            prevent_initial_call=True
        )
        def copy_selected_loci(n_clicks, table_data):
            return self._handle_copy_loci(n_clicks, table_data)

        # MAF and Locus input callbacks
        self._register_filter_callbacks()

    def _handle_download_all_loci(self, selected_samples, selected_subtraits, filter_states, input_loci):
        """Handle download all loci functionality"""
        # MAF 필터 값 추출
        maf_threshold = filter_states.get('maf')
        
        # Collect all sample loci (전체 데이터, trait 필터 없음)
        all_loci = set()
        sample_data_map = {}
        
        # 샘플 데이터 처리
        if selected_samples:
            for sample in selected_samples:
                df_s = self.sample_dfs[sample]
                # trait 필터 제거 - 샘플의 모든 데이터
                df_filt = df_s.copy()
                
                # Apply MAF filter if exists (교집합)
                if maf_threshold is not None:
                    df_filt = df_filt[df_filt['MAF'].astype(float) > maf_threshold]
                
                # Remove duplicates by Locus to ensure unique index
                df_filt = df_filt.drop_duplicates(subset='Locus', keep='first')
                
                sample_loci = set(df_filt['Locus'].tolist())
                all_loci.update(sample_loci)
                
                # Store sample data for lookup
                sample_data_map[sample] = df_filt.set_index('Locus').to_dict('index')
        else:
            # 선택된 샘플이 없으면 빈 sample_data_map 초기화
            for sample in []:  # 빈 리스트
                sample_data_map[sample] = {}
        
        # Input loci 추가 처리
        input_loci_processed = []
        if input_loci:
            for locus in input_loci:
                # chr 접두사 제거해서 원본 형태로 변환
                locus_original = locus.replace('chr', '') if locus.startswith('chr') else locus
                input_loci_processed.append(locus_original)
                all_loci.add(locus_original)
        
        if not all_loci:
            return no_update
        
        # Build table rows for each unique locus (similar to update_locus_table_data)
        rows = []
        for locus in all_loci:
            is_input_locus = locus in input_loci_processed
            
            # 기존 샘플 데이터에서 locus 정보 찾기
            locus_info = None
            if selected_samples:
                for sample in selected_samples:
                    if locus in sample_data_map[sample]:
                        locus_info = sample_data_map[sample][locus]
                        break
            
            # 샘플에 없으면 GWAS에서 가져오기 (input loci의 경우)
            if not locus_info and is_input_locus:
                gwas_data = self.full_df[self.full_df['Locus'] == locus]
                if not gwas_data.empty:
                    locus_info = gwas_data.iloc[0].to_dict()
            
            if not locus_info:
                continue
                
            row_data = {
                'Source': 'inputlocus' if is_input_locus else 'Sample',
                'inputlocus': 'Yes' if is_input_locus else 'No',
                'Variant ID': locus_info.get('Variant ID', ''),
                'Locus': f"chr{locus}",
                'Minor Allele': locus_info.get('Minor Allele', ''),
                'MAF': locus_info.get('MAF', ''),
                'P-value': locus_info.get('P-value', ''),
                'PMID': locus_info.get('PMID', ''),
                'Trait': locus_info.get('Trait', ''),
                'Subtrait': locus_info.get('Subtrait', '')
            }
            
            # Add sample presence and value columns
            if selected_samples:
                for sample in selected_samples:
                    if locus in sample_data_map[sample]:
                        #row_data[f'{sample}_present'] = sample
                        row_data[f'{sample}_value'] = sample_data_map[sample][locus].get('sample_check2', '')
                    else:
                        #row_data[f'{sample}_present'] = ""
                        row_data[f'{sample}_value'] = ""
            
            rows.append(row_data)
        
        if not rows:
            return no_update
        
        # Convert to DataFrame and return as CSV
        df_download = pd.DataFrame(rows)
        return dcc.send_data_frame(df_download.to_csv, "all_sample_loci.csv", index=False)

    def _handle_copy_loci(self, n_clicks, table_data):
        """Handle copy loci functionality"""
        if n_clicks and table_data:
            all_loci = []
            for row in table_data:
                locus = row.get('Locus', '')
                if locus:
                    all_loci.append(locus)
            
            if all_loci:
                copy_text = '\n'.join(all_loci)
                return html.Div([
                    html.Div([
                        html.Span(f"✅ All Loci from the table ({len(all_loci)} items) have been copied to clipboard!", 
                                 style={'font-weight': 'bold', 'color': 'green', 'margin-right': '10px'}),
                        dcc.Clipboard(
                            content=copy_text,
                            title="Copy again",
                            style={"display": "inline-block", "fontSize": 20, "cursor": "pointer"}
                        )
                    ], style={'display': 'flex', 'align-items': 'center', 'margin-bottom': '10px'}),
                    html.Details([
                        html.Summary("Preview copied content", style={'cursor': 'pointer', 'color': 'gray'}),
                        html.Pre(copy_text, style={'background-color': '#f5f5f5', 'padding': '10px', 'border-radius': '4px', 'max-height': '200px', 'overflow-y': 'auto', 'font-family': 'monospace'})
                    ])
                ], style={'padding': '10px', 'border': '1px solid #28a745', 'border-radius': '4px', 'background-color': '#d4edda'})
        elif n_clicks and not table_data:
            return html.Div("No data in table.", style={'color': 'orange', 'font-weight': 'bold'})
        
        return ""

    def _register_filter_callbacks(self):
        """Register MAF and Locus input filter callbacks"""
        # MAF 고려 버튼 클릭 시 입력창 표시
        @self.app.callback(
            [Output('maf-consider-btn', 'style'),
             Output('maf-input-container', 'children')],
            Input('maf-consider-btn', 'n_clicks'),
            prevent_initial_call=True
        )
        def show_maf_input(n_clicks):
            if n_clicks > 0:
                input_component = dcc.Input(
                    id='maf-input',
                    type='number',
                    placeholder='Enter MAF value and press Enter',
                    min=0, max=1, step=0.01,
                    style={'width': '150px', 'margin-right': '10px'}
                )
                return {'display': 'none'}, [input_component]
            return no_update, no_update

        # MAF 입력값 처리 (Enter 또는 blur 이벤트)
        @self.app.callback(
            [Output('maf-input-container', 'children', allow_duplicate=True),
             Output('maf-span-container', 'children'),
             Output('filter-states', 'data', allow_duplicate=True)],
            [Input('maf-input', 'n_submit'),
             Input('maf-input', 'n_blur')],
            [State('maf-input', 'value'),
             State('filter-states', 'data')],
            prevent_initial_call=True
        )
        def process_maf_input(n_submit, n_blur, value, current_states):
            if (n_submit or n_blur) and value is not None:
                # Remove input field
                # Create span
                span_component = html.Span([
                    f"MAF > {value}",
                    html.Button("✕", id="maf-remove-btn", n_clicks=0,
                               style={'margin-left': '5px', 'background': 'transparent', 
                                     'border': 'none', 'color': 'red', 'cursor': 'pointer',
                                     'font-weight': 'bold'})
                ], style={'padding': '5px 10px', 'border': '1px solid #ccc', 
                         'border-radius': '4px', 'background-color': '#f8f9fa',
                         'margin-right': '10px'})
                
                # Update Store state
                updated_states = current_states.copy() if current_states else {}
                updated_states['maf'] = value
                
                return [], [span_component], updated_states
            return no_update, no_update, no_update

        # Locus 입력 버튼 클릭 시 textarea 표시
        @self.app.callback(
            [Output('locus-input-btn', 'style'),
             Output('locus-input-container', 'children')],
            Input('locus-input-btn', 'n_clicks'),
            prevent_initial_call=True
        )
        def show_locus_input(n_clicks):
            if n_clicks > 0:
                input_component = dcc.Textarea(
                    id='locus-textarea',
                    placeholder='Enter loci (e.g. chr3:6340014, chr5:1234567)\nOne per line',
                    style={'width': '300px', 'height': '100px', 'margin-right': '10px'},
                    value=''
                )
                submit_btn = html.Button("Add", id="locus-submit-btn", n_clicks=0,
                                        style={'margin-left': '5px', 'vertical-align': 'top'})
                return {'display': 'none'}, [input_component, submit_btn]
            return no_update, no_update

        # Locus 입력값 처리
        @self.app.callback(
            [Output('locus-input-container', 'children', allow_duplicate=True),
             Output('locus-span-container', 'children'),
             Output('input-loci-store', 'data')],
            Input('locus-submit-btn', 'n_clicks'),
            State('locus-textarea', 'value'),
            prevent_initial_call=True
        )
        def process_locus_input(n_clicks, value):
            if n_clicks > 0 and value:
                # Parse entered loci
                loci_lines = [line.strip() for line in value.split('\n') if line.strip()]
                valid_loci = []
                
                for locus in loci_lines:
                    # Add chr prefix if not present
                    if not locus.startswith('chr'):
                        locus = 'chr' + locus
                    valid_loci.append(locus)
                
                if valid_loci:
                    # Create span
                    span_content = [
                        f"Input Loci ({len(valid_loci)} items)",
                        html.Button("✕", id="locus-remove-btn", n_clicks=0,
                                   style={'margin-left': '5px', 'background': 'transparent', 
                                         'border': 'none', 'color': 'red', 'cursor': 'pointer',
                                         'font-weight': 'bold'})
                    ]
                    span_component = html.Span(span_content, 
                                             style={'padding': '5px 10px', 'border': '1px solid #ccc', 
                                                   'border-radius': '4px', 'background-color': '#e8f4fd',
                                                   'margin-right': '10px'})
                    
                    return [], [span_component], valid_loci
                
            return no_update, no_update, no_update

        # Locus span 제거
        @self.app.callback(
            [Output('locus-span-container', 'children', allow_duplicate=True),
             Output('locus-input-btn', 'style', allow_duplicate=True),
             Output('locus-input-container', 'children', allow_duplicate=True),
             Output('input-loci-store', 'data', allow_duplicate=True)],
            Input('locus-remove-btn', 'n_clicks'),
            prevent_initial_call=True
        )
        def remove_locus_filter(n_clicks):
            if n_clicks > 0:
                # Remove span, reset input field, show button again, reset store
                return [], {'margin-right': '10px', 'margin-left': '20px'}, [], []
            return no_update, no_update, no_update, no_update

    def make_pvalue_scatter_for_samples(self, selected_samples, selected_subtraits=None, maf_threshold=None,
                                        input_loci=None, pvalue_cutoff=None, unique_only=False):
        """
        Enhanced P-value scatter plot with legends and proper chromosome sorting
        
        Args:
            selected_samples: List of sample names
            selected_subtraits: List of subtrait filters
            maf_threshold: MAF filtering threshold
            input_loci: List of input loci
            pvalue_cutoff: P-value cutoff for significance line
            unique_only: Show only sample-unique positions
            
        Returns:
            plotly.graph_objects.Figure: Enhanced scatter plot with legends
        """
        # 데이터 집계
        records = []
        for i, sample in enumerate(selected_samples or []):
            if sample not in self.sample_dfs: 
                continue
            df = self.sample_dfs[sample].copy()
            
            # 컬럼명 맞춤: Start_Position -> Position으로 통일
            if 'Start_Position' in df.columns and 'Position' not in df.columns:
                df['Position'] = df['Start_Position'].astype(int)
            elif 'Position' not in df.columns:
                df['Position'] = 0  # 기본값
            
            if selected_subtraits:
                df = df[df['Subtrait'].isin(selected_subtraits)]
            if maf_threshold is not None:
                df = df[df['MAF'].astype(float) > maf_threshold]
            if input_loci:
                # input_loci는 'chr1:12345' 형태, df의 Locus는 '1:12345' 형태
                normalized_input_loci = []
                for locus in input_loci:
                    if locus.startswith('chr'):
                        normalized_input_loci.append(locus[3:])  # chr 제거
                    else:
                        normalized_input_loci.append(locus)
                df = df[df['Locus'].isin(normalized_input_loci)]
            
            df['sample'] = sample
            df['shape'] = MARKER_SYMBOLS[i % len(MARKER_SYMBOLS)]
            df['color'] = df['Subtrait'].map(lambda sub: SUBTRAIT_COLOR_MAP.get(sub, '#888888'))
            
            # P-value 처리 및 -log10 계산
            df['P-value'] = pd.to_numeric(df['P-value'], errors='coerce')
            valid_pvalue_mask = (df['P-value'] > 0) & pd.notna(df['P-value'])
            df = df[valid_pvalue_mask]
            df['-log10p'] = -np.log10(df['P-value'])
            
            records.append(df)
            
        if not records:
            return go.Figure()
            
        all_df = pd.concat(records, ignore_index=True)

        # unique position 처리
        if unique_only:
            if len(selected_samples) > 1:
                # 여러 샘플: 하나의 샘플에만 있는 position 찾기
                # subtrait filtering이 이미 적용된 상태에서 unique position 계산
                pos_group = all_df.groupby(['Chromosome', 'Position'])['sample'].nunique().reset_index()
                uniq_pos = pos_group[pos_group['sample'] == 1][['Chromosome', 'Position']]
                all_df = all_df.merge(uniq_pos, on=['Chromosome', 'Position'], how='inner')
                
                # subtrait 필터가 적용된 후에도 각 position에 대해 적절한 subtrait만 유지
                if selected_subtraits:
                    all_df = all_df[all_df['Subtrait'].isin(selected_subtraits)]
            else:
                # 단일 샘플: 중복 position 제거 (가장 유의한 P-value만 유지)
                all_df = all_df.loc[all_df.groupby(['Chromosome', 'Position'])['P-value'].idxmin()]
        
        # 마커 심볼 설정
        if not unique_only or len(selected_samples) == 1:
            all_df['shape'] = 'circle'

        # chromosome별 subplot - 고정된 4x3 그리드에 정확한 위치 매핑
        def chr_sort_key(x):
            """염색체를 숫자 순서로 정렬하는 키 함수"""
            if isinstance(x, str):
                if x.startswith('chr'):
                    try:
                        return int(x[3:])  # 'chr' 제거 후 숫자 변환
                    except ValueError:
                        return 999
                else:
                    try:
                        return int(x)  # 숫자 문자열인 경우
                    except ValueError:
                        return 999
            return 999
        
        chromosomes = sorted(all_df['Chromosome'].astype(str).unique(), key=chr_sort_key)
        if len(chromosomes) == 0:
            return go.Figure()
            
        # 고정된 4x3 그리드 설정 (12개 염색체 가정)
        cols = 4
        rows = 3
        
        # subplot titles - 12개 모든 염색체에 대해 설정 (데이터가 없어도 빈 subplot으로 표시)
        subplot_titles = [f'CHR{i}' for i in range(1, 13)]
        
        fig = sp.make_subplots(
            rows=rows, cols=cols,
            subplot_titles=subplot_titles,
            shared_xaxes='columns',
            shared_yaxes='rows',
            vertical_spacing=0.08,
            horizontal_spacing=0.05
        )
        
        # 염색체 번호를 subplot 위치로 매핑하는 딕셔너리 생성
        chr_to_subplot = {}
        for i in range(1, 13):
            row = ((i - 1) // cols) + 1
            col = ((i - 1) % cols) + 1
            chr_to_subplot[str(i)] = (row, col)
            chr_to_subplot[f'chr{i}'] = (row, col)

        # P-value cutoff 수평선
        cutoff_y = None
        if pvalue_cutoff:
            cutoff_y = -np.log10(pvalue_cutoff)

        # subtrait legend 중복 방지를 위한 추적 변수
        subtrait_legend_shown = set()
        
        # 각 염색체별로 플롯 생성 - 정확한 subplot 위치 사용
        for chr_id in chromosomes:
            # 염색체 ID를 정규화하여 subplot 위치 찾기
            chr_normalized = str(chr_id)
            if chr_normalized.startswith('chr'):
                chr_normalized = chr_normalized[3:]
            
            if chr_normalized not in chr_to_subplot and f'chr{chr_normalized}' not in chr_to_subplot:
                continue
                
            # 정확한 subplot 위치 가져오기
            if chr_normalized in chr_to_subplot:
                row, col = chr_to_subplot[chr_normalized]
            else:
                row, col = chr_to_subplot[f'chr{chr_normalized}']
            
            chr_data = all_df[all_df['Chromosome'].astype(str) == str(chr_id)].copy()
            
            if len(chr_data) == 0:
                continue
            
            # Position으로 정렬
            chr_data = chr_data.sort_values('Position')
            
            # subtrait별로 trace 분리하여 legend toggle 기능 구현
            for subtrait in chr_data['Subtrait'].unique():
                subtrait_data = chr_data[chr_data['Subtrait'] == subtrait]
                
                if unique_only and len(selected_samples) > 1:
                    # unique 모드에서는 샘플별로 다른 마커 모양
                    for sample in subtrait_data['sample'].unique():
                        sample_subtrait_data = subtrait_data[subtrait_data['sample'] == sample]
                        
                        # P-value cutoff 적용 시 투명도 조정
                        if pvalue_cutoff and pvalue_cutoff > 0:
                            opacity_values = np.where(sample_subtrait_data['P-value'] <= pvalue_cutoff, 1.0, 0.3)
                        else:
                            opacity_values = [1.0] * len(sample_subtrait_data)
                        
                        fig.add_trace(go.Scattergl(
                            x=sample_subtrait_data['Position'],
                            y=sample_subtrait_data['-log10p'],
                            mode='markers',
                            marker=dict(
                                color=SUBTRAIT_COLOR_MAP.get(subtrait, '#888888'),
                                size=8,
                                symbol=sample_subtrait_data['shape'].iloc[0] if len(sample_subtrait_data) > 0 else 'circle',
                                opacity=opacity_values,
                                line=dict(width=1, color='black')
                            ),
                            text=[f"Pos: {pos:,.0f}<br>-log₁₀(P): {(-np.log10(pv) if pd.notna(pv) and pv > 0 else -1):.1f}<br>Trait: {trait}<br>MAF: {maf}<br>{sample} GT: {gt}<br>Subtrait: {subtr}"
                                   for pos, pv, subtr, trait, maf, gt in zip(sample_subtrait_data['Position'], sample_subtrait_data['P-value'], sample_subtrait_data['Subtrait'], sample_subtrait_data.get('Trait', ['Unknown']*len(sample_subtrait_data)), sample_subtrait_data.get('MAF', ['N/A']*len(sample_subtrait_data)), sample_subtrait_data.get('sample_check2', ['N/A']*len(sample_subtrait_data)))],
                            hovertemplate='%{text}<extra></extra>',
                            hoverinfo='text',
                            name=f"{subtrait}",
                            legendgroup=f"subtrait_{subtrait}",
                            showlegend=(subtrait not in subtrait_legend_shown)
                        ), row=row, col=col)
                        
                        # legend 표시됨을 추적
                        if subtrait not in subtrait_legend_shown:
                            subtrait_legend_shown.add(subtrait)
                else:
                    # 일반 모드에서는 모든 샘플을 원형 마커로
                    # P-value cutoff 적용 시 투명도 조정
                    if pvalue_cutoff and pvalue_cutoff > 0:
                        opacity_values = np.where(subtrait_data['P-value'] <= pvalue_cutoff, 1.0, 0.3)
                    else:
                        opacity_values = [1.0] * len(subtrait_data)
                    
                    fig.add_trace(go.Scattergl(
                        x=subtrait_data['Position'],
                        y=subtrait_data['-log10p'],
                        mode='markers',
                        marker=dict(
                            color=SUBTRAIT_COLOR_MAP.get(subtrait, '#888888'),
                            size=8,
                            symbol='circle',
                            opacity=opacity_values,
                            line=dict(width=1, color='black')
                        ),
                        text=[f"Pos: {pos:,.0f}<br>-log₁₀(P): {(-np.log10(pv) if pd.notna(pv) and pv > 0 else -1):.1f}<br>Trait: {trait}<br>MAF: {maf}<br>{smp} GT: {gt}<br>Subtrait: {subtr}"
                               for smp, pos, pv, subtr, trait, maf, gt in zip(subtrait_data['sample'], subtrait_data['Position'], subtrait_data['P-value'], subtrait_data['Subtrait'], subtrait_data.get('Trait', ['Unknown']*len(subtrait_data)), subtrait_data.get('MAF', ['N/A']*len(subtrait_data)), subtrait_data.get('sample_check2', ['N/A']*len(subtrait_data)))],
                        hovertemplate='%{text}<extra></extra>',
                        hoverinfo='text',
                        name=f"{subtrait}",
                        legendgroup=f"subtrait_{subtrait}",
                        showlegend=(subtrait not in subtrait_legend_shown)
                    ), row=row, col=col)
                    
                    # legend 표시됨을 추적
                    if subtrait not in subtrait_legend_shown:
                        subtrait_legend_shown.add(subtrait)

            # P-value cutoff 수평선 추가
            if cutoff_y is not None:
                fig.add_hline(
                    y=cutoff_y,
                    line_dash="dash",
                    line_color="red",
                    line_width=2,
                    annotation_text=f"-log₁₀(P) = {-np.log10(pvalue_cutoff):.1f}",
                    annotation_position="bottom right",
                    row=row, col=col
                )
            
            # 축 제목 제거 (전체 제목에서 설명)
            fig.update_yaxes(title_text="", row=row, col=col)
            fig.update_xaxes(title_text="", row=row, col=col)

        # 레이아웃 설정 및 범례 추가
        height = int(350 * rows * 1.2)
        fig.update_layout(
            height=height,
            title={
                'text': 'GWAS Variants<br><sub>X-axis: Position (bp) | Y-axis: -log₁₀(P)</sub>',
                'x': 0.5,
                'font': {'size': 16}
            },
            showlegend=True,
            margin=dict(t=100, b=50, l=50, r=50)  # 상단 여백 증가
        )
        
        # unique_only 모드에서 샘플 마커 범례 추가
        if unique_only and len(selected_samples) > 1:
            for i, sample in enumerate(selected_samples):
                marker_symbol = MARKER_SYMBOLS[i % len(MARKER_SYMBOLS)]
                fig.add_trace(
                    go.Scatter(
                        x=[None], y=[None],
                        mode='markers',
                        marker=dict(size=10, color='gray', symbol=marker_symbol),
                        name=f"Sample: {sample}",
                        showlegend=True,
                        legendgroup="samples"
                    )
                )
        
        return fig


# Convenience function for easy usage
def create_gwas_app(data_path=None, 
                   filtered_results_path='./filtered_results_2/matched_results_yesonly_2/',
                   **kwargs):
    """
    Convenience function to create a GWAS module
    
    Args:
        data_path: Path to main GWAS CSV file
        filtered_results_path: Path to filtered results directory
        **kwargs: Additional arguments for app creation
    
    Returns:
        GWASModule: Configured GWAS module instance
    """
    gwas_module = GWASModule(data_path, filtered_results_path)
    gwas_module.create_app(**kwargs)
    return gwas_module


def demo_gwas_module2_usage():
    """
    사용자가 요청한 예시 코드 데모
    """
    # 예시 sample_dfs 생성 (실제로는 get_gwas_data_for_samples에서 가져옴)
    import pandas as pd
    import numpy as np
    
    # sample_dfs = {
    #     'sample1': pd.DataFrame(...),  # 반드시 columns: Chromosome, Position, Subtrait, P-value, MAF, ...
    #     'sample2': pd.DataFrame(...),
    #     ...
    # }
    
    # GWASModule 인스턴스 생성
    gwas_module = GWASModule()
    
    # Bar chart 그대로 (기존 방식)
    # fig_bar = gwas_module.make_figure_for_subtrait(['sample1', 'sample2'])
    
    # Scatter plot - y축 pvalue, unique 모드 OFF
    fig_scatter = gwas_module.make_pvalue_scatter_for_samples(
        selected_samples=['sample1', 'sample2'],
        selected_subtraits=None,
        maf_threshold=0.05,
        input_loci=None,
        pvalue_cutoff=0.01,
        unique_only=False
    )
    
    # Scatter plot - y축 pvalue, unique 모드 ON
    fig_scatter_unique = gwas_module.make_pvalue_scatter_for_samples(
        selected_samples=['sample1', 'sample2'],
        selected_subtraits=None,
        maf_threshold=0.05,
        input_loci=None,
        pvalue_cutoff=0.01,
        unique_only=True
    )
    
    return fig_scatter, fig_scatter_unique 