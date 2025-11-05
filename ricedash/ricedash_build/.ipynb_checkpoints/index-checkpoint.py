import dash
import dash_bootstrap_components as dbc
from dash import html, dcc
from dash.dependencies import Input, Output
import dash_cytoscape as cyto
import pathlib

# 통일된 스타일 import
try:
    from assets.pedi_utils import COMMON_STYLES, create_header, create_footer
except ImportError:
    # 백업 스타일
    COMMON_STYLES = {
        'font_family': '"Nanum Gothic", sans-serif',
        'primary_color': '#4A90E2',
        'background_color': '#FFFFFF'
    }

# 페이지 모듈들 import 추가
import app7_home
import app7_1

# app7_2는 조건부 import
try:
    import app7_2
    APP7_2_AVAILABLE = True
except ImportError as e:
    print(f"⚠️ app7_2 import failed: {e}")
    APP7_2_AVAILABLE = False

# Cytoscape 레이아웃 초기화
cyto.load_extra_layouts()

# app 초기화
app = dash.Dash(__name__, 
                external_stylesheets=[
                    dbc.themes.BOOTSTRAP,
                    'https://use.fontawesome.com/releases/v5.15.4/css/all.css'
                ],
                suppress_callback_exceptions=True)
server = app.server

# index_string 설정
app.index_string = '''
<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>Rice DBreeder</title>
        {%favicon%}
        {%css%}
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>
'''

# 기본 레이아웃 설정
app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content')
], style={
    'backgroundColor': COMMON_STYLES['background_color'],
    'fontFamily': COMMON_STYLES['font_family'],
    'minHeight': '100vh'
})

# 간단한 오류 처리를 포함한 페이지 라우팅
@app.callback(
    Output('page-content', 'children'),
    [Input('url', 'pathname'),
     Input('url', 'search')]
)
def display_page(pathname, search):
    try:
        if pathname == '/app7_1':
            return app7_1.layout
        elif pathname == '/app7_2':
            processed_name = ''
            variety_id = ''
            
            if search:
                try:
                    params = dict(param.split('=') for param in search[1:].split('&'))
                    processed_name = params.get('processed_name', '')
                    variety_id = params.get('id', '')
                    print(f"🎯 index.py URL processing: {params}")
                except Exception as e:
                    print(f"URL 파싱 오류: {e}")
            
            return app7_2.create_layout(processed_name, variety_id)
        elif pathname == '/tutorial':
            md_path = pathlib.Path(__file__).parent / "assets" / "riceDBreeder_tutorial_modern_v4.md"
            with open(md_path, encoding="utf-8") as f:
                md_text = f.read()

            return html.Div([
                create_header(),
                html.Div(className="tut-wrap", children=[
                    html.Div(className="tut-grid", children=[
                        html.Div(className="tut-main", children=[
                            dcc.Markdown(
                                md_text,
                                link_target="_blank",
                                dangerously_allow_html=True,
                            )
                        ]),
                        html.Aside(className="tut-side", children=[
                            html.H4("Key Sections"),
                            html.A("Background & Data Integration", href="#background"),
                            html.A("Overview Preview", href="#preview"),
                            html.A("Using the App", href="#using"),
                            html.Hr(),
                            html.A("– Home Interface", href="#home", style={"marginLeft": "20px", "fontSize": "13px"}),
                            html.A("– Phenotype Search", href="#phenotype-search", style={"marginLeft": "20px", "fontSize": "13px"}),
                            html.A("– Integrated Visual Dashboard", href="#visual-dashboard", style={"marginLeft": "20px", "fontSize": "13px"}),
                            html.Hr(),
                            html.A("– Pedigree Interaction Controls", href="#controls", style={"marginLeft": "20px", "fontSize": "13px"}),
                            html.A("– Add Option Controls", href="#add-option", style={"marginLeft": "20px", "fontSize": "13px"}),
                            html.A("– Phenotype Tab", href="#phenotype", style={"marginLeft": "20px", "fontSize": "13px"}),
                            html.A("– GWAS Tab", href="#gwas", style={"marginLeft": "20px", "fontSize": "13px"}),
                            html.Hr(),
                            html.A("Back to top ↑", href="#top")
                        ]),
                    ])
                ]),
                create_footer()
            ])

        else:
            return app7_home.layout
    except Exception as e:
        print(f"Error loading page {pathname}: {e}")
        import traceback
        traceback.print_exc()
        return html.Div([
            html.H1("Page Load Error"),
            html.P("Unable to load the requested page. Please try again."),
            html.P(f"Error details: {str(e)}"),
            dcc.Link("Return to Home", href="/")
        ])

if __name__ == '__main__':
    app.run_server(debug=False, host='0.0.0.0')
