# pedi_app124.py의 유틸리티 함수들을 여기로 이동
# (normalize_variety_name, fetch_all_resource_ids, normalize_cross_name,
#  fetch_variety_info, try_pedigree_search, fetch_pedigree_data,
#  fetch_pedigree_details, create_cytoscape_elements, default_stylesheet)

# 기존 코드에서 해당 함수들을 그대로 복사하여 이곳에 붙여넣으면 됩니다. 
from dash import dcc

import psycopg2
from psycopg2 import pool
import pandas as pd
import os
from dotenv import load_dotenv
from dash import html

# =============================================================================
# UNIFIED STYLES (통일된 스타일 정의)
# =============================================================================
COMMON_STYLES = {
    'font_family': '"Nanum Gothic", sans-serif',
    'primary_color': '#4A90E2',
    'secondary_color': '#333333',
    'background_color': '#FFFFFF',
    'margin': '20px',
    'padding': '10px'
}

HEADER_STYLE = {
    'backgroundColor': COMMON_STYLES['primary_color'],
    'color': COMMON_STYLES['background_color'],
    'padding': COMMON_STYLES['padding'],
    'textAlign': 'center',
    'marginBottom': COMMON_STYLES['margin']
}

FOOTER_STYLE = {
    'textAlign': 'center',
    'padding': COMMON_STYLES['padding'],
    'marginTop': COMMON_STYLES['margin'],
    'backgroundColor': '#f8f9fa',
    'color': COMMON_STYLES['secondary_color']
}

# =============================================================================
# UNIFIED HEADER AND FOOTER COMPONENTS
# =============================================================================
def create_header():
    """Create unified header component with link to home"""
    return html.Header(
        children=[
            dcc.Link(
                html.H1('RiceDBreeder', className='app-title'),
                href='/',  # 또는 '/home' 등 원하는 라우팅 주소
                style={'textDecoration': 'none', 'color': 'inherit'}
            )
        ],
        className='main-header',
        style=HEADER_STYLE
    )
def create_footer():
    """Create unified footer component"""
    return html.Footer(
        children=[
            html.P('© 2025 RiceDBreeder', className='footer-text'),
            html.P('Gyeongsang National University', className='footer-text')
        ],
        className='main-footer',
        style=FOOTER_STYLE
    )
 # ✅ 한 곳에서만 pool 생성
# Cytoscape 스타일시트
default_stylesheet = [
    {
        'selector': 'node',
        'style': {
            'content': 'data(label)',
            'text-wrap': 'wrap',
            'text-valign': 'center',
            'text-halign': 'center',
            'background-color': '#607D8B',
            'shape': 'ellipse',
            'width': '60px',
            'height': '60px',
            'font-size': '36px',
            'text-margin-y': '-5px',
            'border-color': '#455A64',
            'border-width': '2px',
            'text-outline-color': '#ffffff',
            'text-outline-width': '2px',
            'color': '#000000'
        }
    },
    {
        'selector': '.search-term',
        'style': {
            'background-color': '#EF5350',
            'border-width': '3px',
            'border-color': '#D32F2F',
            'font-weight': 'bold'
        }
    },
    {
        'selector': '.search-node',
        'style': {
            'background-color': '#ff8c00',
            'border-width': '3px',
            'border-color': '#e67e00',
            'font-weight': 'bold'
        }
    },
    {
        'selector': 'edge',
        'style': {
            'curve-style': 'bezier',
            'target-arrow-shape': 'triangle',
            'line-color': '#90A4AE',
            'target-arrow-color': '#78909C',
            'width': '2px',
            'label': 'data(feature)',
            'font-size': '10px',
            'text-rotation': 'autorotate',
            'text-margin-y': '-10px',
            'text-background-color': 'white',
            'text-background-opacity': 1,
            'text-background-padding': '2px'
        }
    }
]

