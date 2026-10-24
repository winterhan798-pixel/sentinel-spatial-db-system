import streamlit as st
import psycopg2
import geopandas as gpd
import folium
from streamlit_folium import st_folium

# 1. 页面基础配置
st.set_page_config(
    page_title="European Sentinel-2 Spatial Query Engine",
    page_icon="🛰️",
    layout="wide"
)

st.title("🛰️ European Sentinel-2 Imagery Spatial Search Engine")
st.markdown("""
Query local **PostGIS** database in real-time to filter and visualize satellite scene coverage footprint across major European cities.
""")

# 2. 数据库连接设置
DB_CONN_PARAMS = {
    "host": "localhost",
    "database": "sentinel_db",
    "user": "postgres",
    "password": "postgres",
    "port": 5432 # 默认为标准的 5432，如在 Colab 转发等环境中请根据实际情况配置
}

# 3. 侧边栏：交互式空间搜索控件
st.sidebar.header("🔍 Spatial Search Parameters")
st.sidebar.write("Adjust the parameters below to filter imagery database in real-time.")

selected_city = st.sidebar.selectbox(
    "Target European City:",
    ["All Cities", "Munich", "Berlin", "Paris", "Rome", "Madrid", "Vienna", "Brussels", "Amsterdam"]
)

max_cloud = st.sidebar.slider(
    "Maximum Allowed Cloud Cover (%)",
    min_value=0.0,
    max_value=15.0,
    value=10.0,
    step=0.5
)

# 4. 构造动态 PostGIS SQL 空间查询
try:
    conn = psycopg2.connect(**DB_CONN_PARAMS)
    
    query = """
        SELECT id, scene_id, acquisition_date, cloud_cover, footprint_geom as geom 
        FROM Scene_Imagery 
        WHERE cloud_cover <= %s
    """
    params = [float(max_cloud)]
    
    if selected_city != "All Cities":
        city_pattern = f"%_{selected_city.strip().upper()[:3]}_%"
        query += " AND scene_id LIKE %s"
        params.append(city_pattern)
        
    # 使用 GeoPandas 提取空间数据并绑定 geom 空间列
    gdf = gpd.read_postgis(query, con=conn, geom_col="geom", params=params)
    conn.close()
    
    # 5. 右侧主屏幕：地图展示与元数据列表
    col1, col2 = st.columns([3, 1])
    
    with col1:
        st.subheader("🗺️ Interactive Imagery Footprints Map")
        if not gdf.empty:
            # 计算中心坐标
            centroid = gdf.geometry.unary_union.centroid
            map_center = [centroid.y, centroid.x]
            
            # 使用清爽且无水印的 Voyager 瓦片地图
            tile_url = 'https://basemaps.cartocdn.com/rastertiles/voyager/{z}/{x}/{y}.png?key=cb1_4tir_1_fe69fb450c80edbef7a3a34e'
            attribution = '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors &copy; <a href="https://carto.com/attributions">CARTO</a>'
            
            m = folium.Map(
                location=map_center, 
                zoom_start=5 if selected_city == "All Cities" else 9, 
                tiles=tile_url,
                attr=attribution,
                control_scale=True
            )
            
            for idx, row in gdf.iterrows():
                is_munich = "_MUN_" in row["scene_id"]
                color = "red" if is_munich else "blue"
                fill_color = "red" if is_munich else "dodgerblue"
                
                popup_html = f"""
                <div style='font-family: Arial, sans-serif; font-size: 11px;'>
                    <b>Scene ID:</b> {row['scene_id']}<br>
                    <b>Acquisition Date:</b> {row['acquisition_date'].strftime('%Y-%m-%d %H:%M')}<br>
                    <b>Cloud Cover:</b> <span style='color: {'red' if row['cloud_cover'] > 10 else 'green'}; font-weight: bold;'>{row['cloud_cover']}%</span>
                </div>
                """
                
                folium.GeoJson(
                    row['geom'],
                    style_function=lambda x, color=color, fill_color=fill_color: {
                        'fillColor': fill_color, 'color': color, 'weight': 2.5 if is_munich else 1.5, 'fillOpacity': 0.15 if is_munich else 0.08
                    },
                    popup=folium.Popup(popup_html, max_width=300)
                ).add_to(m)
            
            # 使用 streamlit_folium 渲染组件
            st_folium(m, width=900, height=550, use_container_width=True)
        else:
            st.warning("⚠️ No satellite scenes match the current filtering parameters. Try increasing the maximum cloud cover or selecting 'All Cities'.")
            
    with col2:
        st.subheader("📊 Statistics & Metadata")
        st.metric(label="Filtered Scenes Count", value=f"{len(gdf)} / 20")
        
        if not gdf.empty:
            st.write("**Metadata List**")
            st.dataframe(
                gdf[['scene_id', 'cloud_cover', 'acquisition_date']].sort_values(by='cloud_cover'),
                height=450,
                hide_index=True
            )
            
except Exception as e:
    st.error(f"❌ Failed to fetch database records: {e}")
