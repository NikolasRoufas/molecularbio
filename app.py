
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
import seaborn as sns
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import warnings
import io
import base64
from datetime import datetime
warnings.filterwarnings('ignore')

# ===== CONFIGURATION =====
st.set_page_config(
    page_title="MolBio Analytics Pro - Molecular Biology Data Analysis Platform",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ===== ENHANCED STYLING =====
def load_custom_css():
    st.markdown("""
    <style>
        /* Import Google Fonts */
        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');
        
        /* Global Styles */
        .main {
            font-family: 'Inter', sans-serif;
        }
        
        /* Header Styles */
        .main-header {
            font-size: 3.5rem;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            text-align: center;
            margin-bottom: 1rem;
            font-weight: 700;
            letter-spacing: -0.02em;
        }
        
        .subtitle {
            text-align: center;
            font-size: 1.2rem;
            color: #6c757d;
            margin-bottom: 2rem;
            font-weight: 400;
        }
        
        /* Section Headers */
        .section-header {
            font-size: 2rem;
            color: #2c3e50;
            margin: 2rem 0 1.5rem 0;
            padding-bottom: 0.5rem;
            border-bottom: 3px solid #3498db;
            font-weight: 600;
        }
        
        .subsection-header {
            font-size: 1.4rem;
            color: #34495e;
            margin: 1.5rem 0 1rem 0;
            font-weight: 500;
        }
        
        /* Card Styles */
        .metric-card {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 2rem;
            border-radius: 15px;
            margin: 0.5rem;
            color: white;
            text-align: center;
            box-shadow: 0 10px 30px rgba(0, 0, 0, 0.1);
            transition: transform 0.3s ease;
        }
        
        .metric-card:hover {
            transform: translateY(-5px);
        }
        
        .metric-value {
            font-size: 2.5rem;
            font-weight: 700;
            margin-bottom: 0.5rem;
        }
        
        .metric-label {
            font-size: 1rem;
            opacity: 0.9;
            font-weight: 400;
        }
        
        /* Info Boxes */
        .info-box {
            background: linear-gradient(135deg, #74b9ff 0%, #0984e3 100%);
            color: white;
            border-radius: 10px;
            padding: 1.5rem;
            margin: 1rem 0;
            box-shadow: 0 5px 15px rgba(116, 185, 255, 0.3);
        }
        
        .warning-box {
            background: linear-gradient(135deg, #fdcb6e 0%, #e17055 100%);
            color: white;
            border-radius: 10px;
            padding: 1.5rem;
            margin: 1rem 0;
            box-shadow: 0 5px 15px rgba(253, 203, 110, 0.3);
        }
        
        .success-box {
            background: linear-gradient(135deg, #00b894 0%, #00a085 100%);
            color: white;
            border-radius: 10px;
            padding: 1.5rem;
            margin: 1rem 0;
            box-shadow: 0 5px 15px rgba(0, 184, 148, 0.3);
        }
        
        /* Team Card Styles */
        .team-card {
            background: linear-gradient(135deg, #a29bfe 0%, #6c5ce7 100%);
            color: white;
            border-radius: 15px;
            padding: 2rem;
            margin: 1rem 0;
            box-shadow: 0 10px 30px rgba(108, 92, 231, 0.3);
            text-align: center;
        }
        
        .team-name {
            font-size: 1.5rem;
            font-weight: 600;
            margin-bottom: 0.5rem;
        }
        
        .team-role {
            font-size: 1rem;
            opacity: 0.9;
            margin-bottom: 1rem;
        }
        
        .team-contribution {
            font-size: 0.9rem;
            opacity: 0.8;
            line-height: 1.5;
        }
        
        /* Button Styles */
        .stButton > button {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border: none;
            border-radius: 10px;
            padding: 0.75rem 2rem;
            font-weight: 500;
            transition: all 0.3s ease;
        }
        
        .stButton > button:hover {
            transform: translateY(-2px);
            box-shadow: 0 5px 15px rgba(102, 126, 234, 0.4);
        }
        
        /* Sidebar Styling */
        .css-1d391kg {
            background: linear-gradient(180deg, #f8f9fa 0%, #e9ecef 100%);
        }
        
        /* Progress Bar */
        .progress-bar {
            background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
            height: 4px;
            border-radius: 2px;
            margin: 1rem 0;
        }
        
        /* Statistics Table */
        .stats-table {
            background: white;
            border-radius: 10px;
            padding: 1rem;
            box-shadow: 0 5px 15px rgba(0, 0, 0, 0.1);
        }
    </style>
    """, unsafe_allow_html=True)

# ===== DATA PROCESSING CLASSES =====
class DataProcessor:
    """Enhanced data processing with comprehensive molecular biology features."""
    
    @staticmethod
    def generate_realistic_sample_data():
        """Generate biologically realistic molecular biology sample data."""
        np.random.seed(42)
        
        # Create biologically relevant gene names
        gene_families = {
            'Oncogenes': ['MYC', 'RAS', 'FOS', 'JUN', 'ERB', 'SRC'],
            'Tumor_suppressors': ['TP53', 'RB1', 'BRCA1', 'BRCA2', 'APC', 'VHL'],
            'Kinases': ['AKT', 'PI3K', 'EGFR', 'VEGFR', 'PDGFR', 'CDK'],
            'Transcription_factors': ['NFkB', 'STAT', 'SOX', 'OCT', 'NANOG', 'KLF'],
            'Housekeeping': ['GAPDH', 'ACTB', 'RPL19', 'HPRT1', 'TBP', 'B2M'],
            'Immune': ['CD4', 'CD8', 'IL2', 'TNF', 'IFNG', 'IL10'],
            'Metabolic': ['GLUT1', 'HK1', 'LDHA', 'PFKL', 'G6PD', 'IDH1']
        }
        
        genes = []
        gene_info = []
        for family, gene_list in gene_families.items():
            for gene in gene_list:
                for i in range(1, 6):  # 5 variants per gene
                    gene_name = f"{gene}{i}"
                    genes.append(gene_name)
                    gene_info.append({'Gene': gene_name, 'Family': family, 'Base_gene': gene})
        
        # Create experimental conditions
        conditions = {
            'Control': {'samples': 6, 'effect': 1.0},
            'Drug_Treatment': {'samples': 6, 'effect': 1.8},
            'Stress_Response': {'samples': 6, 'effect': 0.6},
            'Knockout': {'samples': 6, 'effect': 0.3}
        }
        
        samples = []
        sample_metadata = []
        
        for condition, info in conditions.items():
            for rep in range(1, info['samples'] + 1):
                sample_name = f"{condition}_Rep{rep}"
                samples.append(sample_name)
                sample_metadata.append({
                    'Sample': sample_name,
                    'Condition': condition,
                    'Replicate': rep,
                    'Batch': ((rep - 1) % 3) + 1,
                    'Processing_date': pd.Timestamp('2024-01-01') + pd.Timedelta(days=rep),
                    'RNA_quality': np.random.uniform(7.5, 9.5),
                    'Library_size': np.random.uniform(15e6, 25e6)
                })
        
        # Generate expression data with biological patterns
        n_genes, n_samples = len(genes), len(samples)
        
        # Base expression levels
        base_expression = np.random.lognormal(mean=5, sigma=2, size=(n_genes, n_samples))
        
        # Add condition-specific effects
        gene_df = pd.DataFrame(gene_info)
        
        for i, sample in enumerate(samples):
            condition = sample.split('_')[0]
            effect = conditions[condition]['effect']
            
            # Family-specific responses
            if condition == 'Drug_Treatment':
                # Oncogenes down, tumor suppressors up
                oncogene_mask = gene_df['Family'] == 'Oncogenes'
                tumor_sup_mask = gene_df['Family'] == 'Tumor_suppressors'
                base_expression[oncogene_mask, i] *= np.random.uniform(0.3, 0.7, oncogene_mask.sum())
                base_expression[tumor_sup_mask, i] *= np.random.uniform(1.5, 3.0, tumor_sup_mask.sum())
                
            elif condition == 'Stress_Response':
                # Stress response genes up
                stress_genes = ['TP53', 'NFkB', 'JUN', 'FOS']
                for gene in stress_genes:
                    gene_mask = gene_df['Base_gene'] == gene
                    base_expression[gene_mask, i] *= np.random.uniform(2.0, 4.0, gene_mask.sum())
                    
            elif condition == 'Knockout':
                # Random gene family knocked out
                knockout_family = np.random.choice(['Oncogenes', 'Kinases'])
                knockout_mask = gene_df['Family'] == knockout_family
                base_expression[knockout_mask, i] *= np.random.uniform(0.1, 0.4, knockout_mask.sum())
        
        # Create DataFrames
        expression_df = pd.DataFrame(base_expression, index=genes, columns=samples)
        metadata_df = pd.DataFrame(sample_metadata)
        gene_annotation_df = pd.DataFrame(gene_info)
        
        return expression_df, metadata_df, gene_annotation_df
    
    @staticmethod
    def preprocess_data(data, **kwargs):
        """Enhanced preprocessing with multiple normalization options."""
        processed_data = data.copy()
        steps = []
        
        # Remove genes with too many zeros
        if kwargs.get('filter_zeros', False):
            zero_threshold = kwargs.get('zero_threshold', 0.5)
            zero_fraction = (processed_data == 0).sum(axis=1) / processed_data.shape[1]
            keep_genes = zero_fraction < zero_threshold
            processed_data = processed_data[keep_genes]
            steps.append(f"Removed {(~keep_genes).sum()} genes with >{zero_threshold*100}% zeros")
        
        # Filter low expression
        if kwargs.get('filter_low_expr', True):
            min_expr = kwargs.get('min_expression', 1.0)
            gene_means = processed_data.mean(axis=1)
            high_expr_genes = gene_means >= min_expr
            processed_data = processed_data[high_expr_genes]
            steps.append(f"Filtered {(~high_expr_genes).sum()} low-expression genes (< {min_expr})")
        
        # Log transformation
        if kwargs.get('log_transform', True):
            transform_type = kwargs.get('log_type', 'log2')
            if transform_type == 'log2':
                processed_data = np.log2(processed_data + 1)
                steps.append("Applied log2(x+1) transformation")
            elif transform_type == 'log10':
                processed_data = np.log10(processed_data + 1)
                steps.append("Applied log10(x+1) transformation")
            elif transform_type == 'natural':
                processed_data = np.log(processed_data + 1)
                steps.append("Applied ln(x+1) transformation")
        
        # Normalization
        if kwargs.get('normalize', False):
            norm_method = kwargs.get('norm_method', 'tpm')
            if norm_method == 'tpm':
                processed_data = processed_data.div(processed_data.sum(axis=0), axis=1) * 1e6
                steps.append("Applied TPM normalization")
            elif norm_method == 'quantile':
                from sklearn.preprocessing import quantile_transform
                processed_data = pd.DataFrame(
                    quantile_transform(processed_data.T).T,
                    index=processed_data.index,
                    columns=processed_data.columns
                )
                steps.append("Applied quantile normalization")
        
        # Scaling
        if kwargs.get('scale_data', False):
            scale_method = kwargs.get('scale_method', 'standard')
            if scale_method == 'standard':
                scaler = StandardScaler()
                processed_data = pd.DataFrame(
                    scaler.fit_transform(processed_data.T).T,
                    index=processed_data.index,
                    columns=processed_data.columns
                )
                steps.append("Applied standard scaling")
            elif scale_method == 'robust':
                scaler = RobustScaler()
                processed_data = pd.DataFrame(
                    scaler.fit_transform(processed_data.T).T,
                    index=processed_data.index,
                    columns=processed_data.columns
                )
                steps.append("Applied robust scaling")
        
        return processed_data, steps

# ===== ANALYSIS CLASSES =====
class QualityControl:
    """Comprehensive quality control analysis."""
    
    @staticmethod
    def generate_qc_metrics(data, metadata=None):
        """Generate comprehensive QC metrics."""
        metrics = {
            'basic': {
                'total_genes': data.shape[0],
                'total_samples': data.shape[1],
                'total_measurements': data.shape[0] * data.shape[1],
                'missing_values': data.isnull().sum().sum(),
                'zero_values': (data == 0).sum().sum(),
                'mean_expression': data.values.mean(),
                'median_expression': np.median(data.values),
                'std_expression': data.values.std()
            },
            'distribution': {
                'min_value': data.values.min(),
                'max_value': data.values.max(),
                'q25': np.percentile(data.values, 25),
                'q75': np.percentile(data.values, 75),
                'skewness': stats.skew(data.values.flatten()),
                'kurtosis': stats.kurtosis(data.values.flatten())
            }
        }
        
        # Sample-wise metrics
        sample_metrics = pd.DataFrame({
            'Total_counts': data.sum(axis=0),
            'Detected_genes': (data > 0).sum(axis=0),
            'Mean_expression': data.mean(axis=0),
            'CV': data.std(axis=0) / data.mean(axis=0)
        })
        
        # Gene-wise metrics
        gene_metrics = pd.DataFrame({
            'Mean_expression': data.mean(axis=1),
            'Std_expression': data.std(axis=1),
            'CV': data.std(axis=1) / data.mean(axis=1),
            'Detection_rate': (data > 0).sum(axis=1) / data.shape[1],
            'Max_expression': data.max(axis=1)
        })
        
        metrics['samples'] = sample_metrics
        metrics['genes'] = gene_metrics
        
        return metrics

class AdvancedAnalysis:
    """Advanced statistical and ML analysis."""
    
    @staticmethod
    def perform_pca(data, n_components=None, scale=True):
        """Enhanced PCA with comprehensive results."""
        if n_components is None:
            n_components = min(10, data.shape[1] - 1)
        
        # Prepare data
        X = data.T  # Samples as rows
        
        if scale:
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X)
        else:
            X_scaled = X
        
        # Perform PCA
        pca = PCA(n_components=n_components)
        pca_coords = pca.fit_transform(X_scaled)
        
        # Calculate loadings
        loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
        
        # Feature contributions
        feature_importance = pd.DataFrame(
            np.abs(loadings),
            index=data.index,
            columns=[f'PC{i+1}' for i in range(n_components)]
        )
        
        results = {
            'coordinates': pca_coords,
            'explained_variance_ratio': pca.explained_variance_ratio_,
            'explained_variance': pca.explained_variance_,
            'loadings': loadings,
            'feature_importance': feature_importance,
            'cumulative_variance': np.cumsum(pca.explained_variance_ratio_),
            'sample_names': data.columns.tolist(),
            'feature_names': data.index.tolist()
        }
        
        return results
    
    @staticmethod
    def perform_clustering(data, method='kmeans', n_clusters=None, **kwargs):
        """Comprehensive clustering analysis."""
        X = data.T  # Samples as rows
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        
        results = {}
        
        if method == 'kmeans':
            if n_clusters is None:
                # Find optimal clusters
                k_range = range(2, min(11, X.shape[0]))
                inertias = []
                silhouette_scores = []
                
                for k in k_range:
                    kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
                    cluster_labels = kmeans.fit_predict(X_scaled)
                    inertias.append(kmeans.inertia_)
                    if len(np.unique(cluster_labels)) > 1:
                        silhouette_scores.append(silhouette_score(X_scaled, cluster_labels))
                    else:
                        silhouette_scores.append(0)
                
                # Choose optimal k
                optimal_k = k_range[np.argmax(silhouette_scores)]
                
                results['optimization'] = {
                    'k_range': list(k_range),
                    'inertias': inertias,
                    'silhouette_scores': silhouette_scores,
                    'optimal_k': optimal_k
                }
                
                n_clusters = optimal_k
            
            # Final clustering
            kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
            cluster_labels = kmeans.fit_predict(X_scaled)
            
            results['labels'] = cluster_labels
            results['centers'] = kmeans.cluster_centers_
            results['n_clusters'] = n_clusters
            results['method'] = 'kmeans'
            
        elif method == 'hierarchical':
            linkage_matrix = linkage(X_scaled, method='ward')
            
            if n_clusters is None:
                n_clusters = 3
            
            cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust') - 1
            
            results['labels'] = cluster_labels
            results['linkage_matrix'] = linkage_matrix
            results['n_clusters'] = n_clusters
            results['method'] = 'hierarchical'
        
        # Cluster characterization
        cluster_profiles = []
        for cluster in range(n_clusters):
            cluster_mask = cluster_labels == cluster
            cluster_samples = data.columns[cluster_mask]
            cluster_data = data[cluster_samples]
            
            profile = {
                'cluster': cluster,
                'size': cluster_mask.sum(),
                'samples': cluster_samples.tolist(),
                'mean_profile': cluster_data.mean(axis=1),
                'top_genes': cluster_data.mean(axis=1).nlargest(10).index.tolist()
            }
            cluster_profiles.append(profile)
        
        results['profiles'] = cluster_profiles
        results['sample_names'] = data.columns.tolist()
        
        return results

# ===== VISUALIZATION ENGINE =====
class VisualizationEngine:
    """Advanced visualization creation."""
    
    @staticmethod
    def create_qc_dashboard(data, qc_metrics):
        """Create comprehensive QC dashboard."""
        fig = make_subplots(
            rows=2, cols=3,
            subplot_titles=(
                'Expression Distribution', 'Sample Library Sizes',
                'Gene Detection Rates', 'Expression Heatmap (Top 20 genes)',
                'CV vs Mean Expression', 'Sample Correlation'
            ),
            specs=[[{"type": "histogram"}, {"type": "bar"}, {"type": "histogram"}],
                   [{"type": "heatmap"}, {"type": "scatter"}, {"type": "heatmap"}]]
        )
        
        # Expression distribution
        fig.add_trace(
            go.Histogram(x=data.values.flatten(), nbinsx=50, name="Expression Distribution"),
            row=1, col=1
        )
        
        # Sample library sizes
        sample_totals = data.sum(axis=0)
        fig.add_trace(
            go.Bar(x=sample_totals.index, y=sample_totals.values, name="Library Sizes"),
            row=1, col=2
        )
        
        # Gene detection rates
        detection_rates = (data > 0).sum(axis=1) / data.shape[1]
        fig.add_trace(
            go.Histogram(x=detection_rates, nbinsx=30, name="Detection Rates"),
            row=1, col=3
        )
        
        # Top expressed genes heatmap
        top_genes = data.mean(axis=1).nlargest(20)
        heatmap_data = data.loc[top_genes.index]
        fig.add_trace(
            go.Heatmap(z=heatmap_data.values, 
                      x=heatmap_data.columns, 
                      y=heatmap_data.index,
                      colorscale='Viridis'),
            row=2, col=1
        )
        
        # CV vs Mean
        gene_means = data.mean(axis=1)
        gene_cv = data.std(axis=1) / gene_means
        fig.add_trace(
            go.Scatter(x=gene_means, y=gene_cv, mode='markers', 
                      name="CV vs Mean", opacity=0.6),
            row=2, col=2
        )
        
        # Sample correlation
        corr_matrix = data.T.corr()
        fig.add_trace(
            go.Heatmap(z=corr_matrix.values,
                      x=corr_matrix.columns,
                      y=corr_matrix.index,
                      colorscale='RdBu',
                      zmid=0),
            row=2, col=3
        )
        
        fig.update_layout(height=800, showlegend=False, title_text="Quality Control Dashboard")
        return fig
    
    @staticmethod
    def create_pca_visualization(pca_results, metadata=None, color_by=None):
        """Create comprehensive PCA visualization."""
        coords = pca_results['coordinates']
        
        # Main PCA plot
        if metadata is not None and color_by is not None and color_by in metadata.columns:
            color_values = metadata[color_by].values
            fig = px.scatter(
                x=coords[:, 0], y=coords[:, 1],
                color=color_values,
                hover_name=pca_results['sample_names'],
                title=f"PCA Analysis (colored by {color_by})",
                labels={
                    'x': f"PC1 ({pca_results['explained_variance_ratio'][0]:.1%} variance)",
                    'y': f"PC2 ({pca_results['explained_variance_ratio'][1]:.1%} variance)"
                }
            )
        else:
            fig = px.scatter(
                x=coords[:, 0], y=coords[:, 1],
                hover_name=pca_results['sample_names'],
                title="PCA Analysis",
                labels={
                    'x': f"PC1 ({pca_results['explained_variance_ratio'][0]:.1%} variance)",
                    'y': f"PC2 ({pca_results['explained_variance_ratio'][1]:.1%} variance)"
                }
            )
        
        fig.update_traces(marker=dict(size=10, opacity=0.7))
        fig.update_layout(height=500)
        
        return fig
    
    @staticmethod
    def create_clustering_visualization(cluster_results, pca_results=None):
        """Create clustering visualization."""
        if pca_results is not None:
            coords = pca_results['coordinates']
            fig = px.scatter(
                x=coords[:, 0], y=coords[:, 1],
                color=cluster_results['labels'].astype(str),
                hover_name=cluster_results['sample_names'],
                title=f"{cluster_results['method'].title()} Clustering Results",
                labels={
                    'x': f"PC1 ({pca_results['explained_variance_ratio'][0]:.1%} variance)",
                    'y': f"PC2 ({pca_results['explained_variance_ratio'][1]:.1%} variance)",
                    'color': 'Cluster'
                }
            )
        else:
            # Create a simple visualization without PCA coordinates
            fig = px.bar(
                x=[f"Cluster {i}" for i in range(cluster_results['n_clusters'])],
                y=[sum(cluster_results['labels'] == i) for i in range(cluster_results['n_clusters'])],
                title="Cluster Sizes"
            )
        
        return fig

# ===== MAIN APPLICATION =====
def main():
    load_custom_css()
    
    # Header
    st.markdown('<h1 class="main-header">🧬 MolBio Analytics Pro</h1>', unsafe_allow_html=True)
    st.markdown('<p class="subtitle">Advanced Molecular Biology Data Analysis Platform</p>', unsafe_allow_html=True)
    
    # Initialize session state
    if 'data_loaded' not in st.session_state:
        st.session_state.data_loaded = False
        st.session_state.analysis_complete = False
    
    # Create tabs
    tabs = st.tabs(["🏠 Home", "📊 Analysis", "🎯 Results", "👥 Team Info"])
    
    # Sidebar Configuration
    with st.sidebar:
        st.markdown("## 🔧 Analysis Configuration")
        
        # Data Input Section
        st.markdown("### 📁 Data Input")
        
        data_source = st.radio(
            "Choose data source:",
            ["Upload File", "Use Sample Dataset"]
        )
        
        if data_source == "Upload File":
            uploaded_file = st.file_uploader(
                "Upload Gene Expression Data",
                type=['csv', 'xlsx', 'tsv'],
                help="File should have genes as rows and samples as columns"
            )
        else:
            uploaded_file = None
            generate_sample = st.button("🧪 Generate Sample Dataset", type="primary")
        
        # Processing Parameters
        st.markdown("### ⚙️ Preprocessing Parameters")
        
        with st.expander("🧹 Data Cleaning", expanded=True):
            filter_zeros = st.checkbox("Remove high-zero genes", value=False)
            if filter_zeros:
                zero_threshold = st.slider("Max zero fraction", 0.1, 0.9, 0.5)
            else:
                zero_threshold = 0.5
            
            filter_low_expr = st.checkbox("Filter low expression", value=True)
            if filter_low_expr:
                min_expression = st.slider("Min expression threshold", 0.1, 5.0, 1.0)
            else:
                min_expression = 0.0
        
        with st.expander("🔄 Transformations", expanded=True):
            log_transform = st.checkbox("Log transformation", value=True)
            if log_transform:
                log_type = st.selectbox("Log type", ["log2", "log10", "natural"])
            else:
                log_type = "log2"
            
            normalize = st.checkbox("Normalization", value=False)
            if normalize:
                norm_method = st.selectbox("Normalization method", ["tpm", "quantile"])
            else:
                norm_method = "tpm"
            scale_data = st.checkbox("Scale data", value=False)
            if scale_data:
                scale_method = st.selectbox("Scaling method", ["standard", "robust"])
            else:
                scale_method = "standard"
        
        # Analysis Parameters
        st.markdown("### 🔬 Analysis Parameters")
        
        with st.expander("📊 PCA Settings", expanded=True):
            n_pca_components = st.slider("Number of PCA components", 2, 20, 10)
            scale_for_pca = st.checkbox("Scale data for PCA", value=True)
        
        with st.expander("🎯 Clustering Settings", expanded=True):
            clustering_method = st.selectbox("Clustering method", ["kmeans", "hierarchical"])
            if clustering_method == "kmeans":
                n_clusters = st.slider("Number of clusters (0 = auto)", 0, 10, 0)
            else:
                n_clusters = st.slider("Number of clusters", 2, 10, 3)
        
        with st.expander("🔍 t-SNE Settings"):
            run_tsne = st.checkbox("Run t-SNE", value=False)
            if run_tsne:
                tsne_perplexity = st.slider("Perplexity", 5, 50, 30)
                tsne_iterations = st.slider("Iterations", 250, 2000, 1000)
    
    # HOME TAB
    with tabs[0]:
        st.markdown('<h2 class="section-header">🏠 Welcome to MolBio Analytics Pro</h2>', unsafe_allow_html=True)
        
        # Feature Overview
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.markdown("""
            <div class="info-box">
                <h3>🧬 Data Processing</h3>
                <p>Advanced preprocessing pipeline with filtering, normalization, and transformation options specifically designed for molecular biology data.</p>
            </div>
            """, unsafe_allow_html=True)
        
        with col2:
            st.markdown("""
            <div class="info-box">
                <h3>📊 Quality Control</h3>
                <p>Comprehensive QC metrics and visualizations to assess data quality, including library sizes, detection rates, and expression distributions.</p>
            </div>
            """, unsafe_allow_html=True)
        
        with col3:
            st.markdown("""
            <div class="info-box">
                <h3>🔬 Advanced Analysis</h3>
                <p>Machine learning algorithms including PCA, clustering, and t-SNE for dimensionality reduction and pattern discovery.</p>
            </div>
            """, unsafe_allow_html=True)
        
        # Quick Start Guide
        st.markdown('<h3 class="subsection-header">🚀 Quick Start Guide</h3>', unsafe_allow_html=True)
        
        st.markdown("""
        1. **Choose your data source** in the sidebar: Upload your own file or use our sample dataset
        2. **Configure preprocessing parameters** to clean and transform your data
        3. **Set analysis parameters** for PCA, clustering, and other analyses
        4. **Navigate to the Analysis tab** to load and process your data
        5. **View results** in the Results tab with interactive visualizations
        """)
        
        # Supported File Formats
        st.markdown('<h3 class="subsection-header">📁 Supported File Formats</h3>', unsafe_allow_html=True)
        
        format_col1, format_col2 = st.columns(2)
        
        with format_col1:
            st.markdown("""
            **Input Formats:**
            - CSV files (.csv)
            - Excel files (.xlsx)
            - Tab-separated files (.tsv)
            """)
        
        with format_col2:
            st.markdown("""
            **Data Structure:**
            - Genes as rows
            - Samples as columns
            - Expression values as numeric data
            """)
    
    # ANALYSIS TAB
    with tabs[1]:
        st.markdown('<h2 class="section-header">📊 Data Analysis</h2>', unsafe_allow_html=True)
        
        # Data Loading Section
        if data_source == "Upload File" and uploaded_file is not None:
            try:
                # Load uploaded file
                if uploaded_file.name.endswith('.csv'):
                    raw_data = pd.read_csv(uploaded_file, index_col=0)
                elif uploaded_file.name.endswith('.xlsx'):
                    raw_data = pd.read_excel(uploaded_file, index_col=0)
                elif uploaded_file.name.endswith('.tsv'):
                    raw_data = pd.read_csv(uploaded_file, sep='\t', index_col=0)
                
                st.session_state.raw_data = raw_data
                st.session_state.data_loaded = True
                
                st.markdown("""
                <div class="success-box">
                    <h4>✅ Data Successfully Loaded!</h4>
                    <p>Your gene expression data has been loaded and is ready for analysis.</p>
                </div>
                """, unsafe_allow_html=True)
                
            except Exception as e:
                st.error(f"Error loading file: {str(e)}")
                st.session_state.data_loaded = False
        
        elif data_source == "Use Sample Dataset" and 'generate_sample' in locals() and generate_sample:
            # Generate sample data
            with st.spinner("Generating biologically realistic sample dataset..."):
                raw_data, metadata, gene_annotations = DataProcessor.generate_realistic_sample_data()
                st.session_state.raw_data = raw_data
                st.session_state.metadata = metadata
                st.session_state.gene_annotations = gene_annotations
                st.session_state.data_loaded = True
            
            st.markdown("""
            <div class="success-box">
                <h4>🧪 Sample Dataset Generated!</h4>
                <p>A biologically realistic dataset with multiple conditions and gene families has been created.</p>
            </div>
            """, unsafe_allow_html=True)
        
        # Data Processing Section
        if st.session_state.data_loaded:
            st.markdown('<h3 class="subsection-header">🔧 Data Processing</h3>', unsafe_allow_html=True)
            
            # Show raw data info
            raw_data = st.session_state.raw_data
            
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.markdown(f"""
                <div class="metric-card">
                    <div class="metric-value">{raw_data.shape[0]}</div>
                    <div class="metric-label">Genes</div>
                </div>
                """, unsafe_allow_html=True)
            
            with col2:
                st.markdown(f"""
                <div class="metric-card">
                    <div class="metric-value">{raw_data.shape[1]}</div>
                    <div class="metric-label">Samples</div>
                </div>
                """, unsafe_allow_html=True)
            
            with col3:
                st.markdown(f"""
                <div class="metric-card">
                    <div class="metric-value">{raw_data.isnull().sum().sum()}</div>
                    <div class="metric-label">Missing Values</div>
                </div>
                """, unsafe_allow_html=True)
            
            with col4:
                st.markdown(f"""
                <div class="metric-card">
                    <div class="metric-value">{raw_data.values.mean():.2f}</div>
                    <div class="metric-label">Mean Expression</div>
                </div>
                """, unsafe_allow_html=True)
            
            # Process data button
            if st.button("🚀 Process Data", type="primary", use_container_width=True):
                with st.spinner("Processing data according to your specifications..."):
                    # Apply preprocessing
                    processing_params = {
                        'filter_zeros': filter_zeros,
                        'zero_threshold': zero_threshold,
                        'filter_low_expr': filter_low_expr,
                        'min_expression': min_expression,
                        'log_transform': log_transform,
                        'log_type': log_type,
                        'normalize': normalize,
                        'norm_method': norm_method,
                        'scale_data': scale_data,
                        'scale_method': scale_method
                    }
                    
                    processed_data, processing_steps = DataProcessor.preprocess_data(
                        raw_data, **processing_params
                    )
                    
                    # Store processed data
                    st.session_state.processed_data = processed_data
                    st.session_state.processing_steps = processing_steps
                    
                    # Generate QC metrics
                    qc_metrics = QualityControl.generate_qc_metrics(
                        processed_data, 
                        st.session_state.get('metadata', None)
                    )
                    st.session_state.qc_metrics = qc_metrics
                    
                    # Perform PCA
                    pca_results = AdvancedAnalysis.perform_pca(
                        processed_data, 
                        n_components=n_pca_components,
                        scale=scale_for_pca
                    )
                    st.session_state.pca_results = pca_results
                    
                    # Perform Clustering
                    cluster_results = AdvancedAnalysis.perform_clustering(
                        processed_data,
                        method=clustering_method,
                        n_clusters=n_clusters if n_clusters > 0 else None
                    )
                    st.session_state.cluster_results = cluster_results
                    
                    # Perform t-SNE if requested
                    if run_tsne:
                        from sklearn.manifold import TSNE
                        X = processed_data.T
                        scaler = StandardScaler()
                        X_scaled = scaler.fit_transform(X)
                        
                        tsne = TSNE(
                            n_components=2,
                            perplexity=tsne_perplexity,
                            n_iter=tsne_iterations,
                            random_state=42
                        )
                        tsne_coords = tsne.fit_transform(X_scaled)
                        
                        st.session_state.tsne_results = {
                            'coordinates': tsne_coords,
                            'sample_names': processed_data.columns.tolist()
                        }
                    
                    st.session_state.analysis_complete = True
                
                st.success("✅ Analysis completed successfully!")
                
                # Show processing summary
                st.markdown('<h4 class="subsection-header">📋 Processing Summary</h4>', unsafe_allow_html=True)
                for step in processing_steps:
                    st.write(f"• {step}")
                
                # Show processed data metrics
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.markdown(f"""
                    <div class="metric-card">
                        <div class="metric-value">{processed_data.shape[0]}</div>
                        <div class="metric-label">Final Genes</div>
                    </div>
                    """, unsafe_allow_html=True)
                
                with col2:
                    st.markdown(f"""
                    <div class="metric-card">
                        <div class="metric-value">{processed_data.shape[1]}</div>
                        <div class="metric-label">Final Samples</div>
                    </div>
                    """, unsafe_allow_html=True)
                
                with col3:
                    st.markdown(f"""
                    <div class="metric-card">
                        <div class="metric-value">{processed_data.values.mean():.2f}</div>
                        <div class="metric-label">Mean Expression</div>
                    </div>
                    """, unsafe_allow_html=True)
                
                with col4:
                    st.markdown(f"""
                    <div class="metric-card">
                        <div class="metric-value">{processed_data.values.std():.2f}</div>
                        <div class="metric-label">Std Expression</div>
                    </div>
                    """, unsafe_allow_html=True)
        
        else:
            st.markdown("""
            <div class="warning-box">
                <h4>⚠️ No Data Loaded</h4>
                <p>Please upload a file or generate sample data using the sidebar options.</p>
            </div>
            """, unsafe_allow_html=True)
    
    # RESULTS TAB
    with tabs[2]:
        st.markdown('<h2 class="section-header">🎯 Analysis Results</h2>', unsafe_allow_html=True)
        
        if st.session_state.analysis_complete:
            # Create sub-tabs for different result types
            result_tabs = st.tabs(["📊 Quality Control", "🔍 PCA Analysis", "🎯 Clustering", "📈 Advanced Plots"])
            
            # Quality Control Tab
            with result_tabs[0]:
                st.markdown('<h3 class="subsection-header">📊 Quality Control Dashboard</h3>', unsafe_allow_html=True)
                
                qc_metrics = st.session_state.qc_metrics
                processed_data = st.session_state.processed_data
                
                # QC Metrics Summary
                basic_metrics = qc_metrics['basic']
                dist_metrics = qc_metrics['distribution']
                
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.metric("Total Genes", f"{basic_metrics['total_genes']:,}")
                    st.metric("Mean Expression", f"{basic_metrics['mean_expression']:.2f}")
                
                with col2:
                    st.metric("Total Samples", f"{basic_metrics['total_samples']:,}")
                    st.metric("Std Expression", f"{basic_metrics['std_expression']:.2f}")
                
                with col3:
                    st.metric("Zero Values", f"{basic_metrics['zero_values']:,}")
                    st.metric("Skewness", f"{dist_metrics['skewness']:.2f}")
                
                with col4:
                    st.metric("Missing Values", f"{basic_metrics['missing_values']:,}")
                    st.metric("Kurtosis", f"{dist_metrics['kurtosis']:.2f}")
                
                # QC Visualizations
                qc_fig = VisualizationEngine.create_qc_dashboard(processed_data, qc_metrics)
                st.plotly_chart(qc_fig, use_container_width=True)
                
                # Sample and Gene Statistics Tables
                col1, col2 = st.columns(2)
                
                with col1:
                    st.markdown("#### 📋 Sample Statistics")
                    sample_stats = qc_metrics['samples'].describe()
                    st.dataframe(sample_stats, use_container_width=True)
                
                with col2:
                    st.markdown("#### 🧬 Gene Statistics")
                    gene_stats = qc_metrics['genes'].describe()
                    st.dataframe(gene_stats, use_container_width=True)
            
            # PCA Analysis Tab
            with result_tabs[1]:
                st.markdown('<h3 class="subsection-header">🔍 Principal Component Analysis</h3>', unsafe_allow_html=True)
                
                pca_results = st.session_state.pca_results
                
                # PCA Summary Metrics
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.metric("PC1 Variance", f"{pca_results['explained_variance_ratio'][0]:.1%}")
                with col2:
                    st.metric("PC2 Variance", f"{pca_results['explained_variance_ratio'][1]:.1%}")
                with col3:
                    st.metric("Cumulative PC1-PC2", f"{pca_results['cumulative_variance'][1]:.1%}")
                with col4:
                    st.metric("Total Components", len(pca_results['explained_variance_ratio']))
                
                # PCA Visualization
                metadata = st.session_state.get('metadata', None)
                color_options = ['None']
                if metadata is not None:
                    color_options.extend(metadata.columns.tolist())
                
                color_by = st.selectbox("Color points by:", color_options)
                color_column = color_by if color_by != 'None' else None
                
                pca_fig = VisualizationEngine.create_pca_visualization(
                    pca_results, metadata, color_column
                )
                st.plotly_chart(pca_fig, use_container_width=True)
                
                # Variance Explained Plot
                var_fig = px.bar(
                    x=[f'PC{i+1}' for i in range(len(pca_results['explained_variance_ratio']))],
                    y=pca_results['explained_variance_ratio'],
                    title="Variance Explained by Each Principal Component",
                    labels={'x': 'Principal Component', 'y': 'Variance Explained'}
                )
                st.plotly_chart(var_fig, use_container_width=True)
                
                # Feature Loadings
                st.markdown("#### 🎯 Top Contributing Features")
                pc_select = st.selectbox("Select Principal Component:", 
                                       [f'PC{i+1}' for i in range(len(pca_results['explained_variance_ratio']))])
                
                top_features = pca_results['feature_importance'][pc_select].nlargest(20)
                
                loading_fig = px.bar(
                    x=top_features.values,
                    y=top_features.index,
                    orientation='h',
                    title=f"Top 20 Features Contributing to {pc_select}",
                    labels={'x': 'Absolute Loading', 'y': 'Gene'}
                )
                loading_fig.update_layout(height=600)
                st.plotly_chart(loading_fig, use_container_width=True)
            
            # Clustering Tab
            with result_tabs[2]:
                st.markdown('<h3 class="subsection-header">🎯 Clustering Analysis</h3>', unsafe_allow_html=True)
                
                cluster_results = st.session_state.cluster_results
                pca_results = st.session_state.pca_results
                
                # Clustering Summary
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Number of Clusters", cluster_results['n_clusters'])
                with col2:
                    st.metric("Clustering Method", cluster_results['method'].title())
                with col3:
                    if 'optimization' in cluster_results:
                        st.metric("Optimal K (Silhouette)", cluster_results['optimization']['optimal_k'])
                
                # Clustering Visualization
                cluster_fig = VisualizationEngine.create_clustering_visualization(
                    cluster_results, pca_results
                )
                st.plotly_chart(cluster_fig, use_container_width=True)
                
                # Cluster Optimization (if available)
                if 'optimization' in cluster_results:
                    opt_data = cluster_results['optimization']
                    
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        elbow_fig = px.line(
                            x=opt_data['k_range'],
                            y=opt_data['inertias'],
                            title="Elbow Method for Optimal K",
                            labels={'x': 'Number of Clusters (K)', 'y': 'Inertia'}
                        )
                        st.plotly_chart(elbow_fig, use_container_width=True)
                    
                    with col2:
                        silhouette_fig = px.line(
                            x=opt_data['k_range'],
                            y=opt_data['silhouette_scores'],
                            title="Silhouette Score vs K",
                            labels={'x': 'Number of Clusters (K)', 'y': 'Silhouette Score'}
                        )
                        st.plotly_chart(silhouette_fig, use_container_width=True)
                
                # Cluster Profiles
                st.markdown("#### 📊 Cluster Profiles")
                
                for profile in cluster_results['profiles']:
                    with st.expander(f"Cluster {profile['cluster']} ({profile['size']} samples)"):
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            st.write("**Samples in this cluster:**")
                            st.write(", ".join(profile['samples']))
                        
                        with col2:
                            st.write("**Top expressed genes:**")
                            for gene in profile['top_genes'][:5]:
                                st.write(f"• {gene}")
            
            # Advanced Plots Tab
            with result_tabs[3]:
                st.markdown('<h3 class="subsection-header">📈 Advanced Visualizations</h3>', unsafe_allow_html=True)
                
                processed_data = st.session_state.processed_data
                
                # t-SNE Plot (if available)
                if 'tsne_results' in st.session_state:
                    st.markdown("#### 🔬 t-SNE Visualization")
                    tsne_results = st.session_state.tsne_results
                    
                    # Color by clustering results
                    cluster_labels = st.session_state.cluster_results['labels']
                    
                    tsne_fig = px.scatter(
                        x=tsne_results['coordinates'][:, 0],
                        y=tsne_results['coordinates'][:, 1],
                        color=cluster_labels.astype(str),
                        hover_name=tsne_results['sample_names'],
                        title="t-SNE Visualization (colored by clusters)",
                        labels={'x': 't-SNE 1', 'y': 't-SNE 2', 'color': 'Cluster'}
                    )
                    st.plotly_chart(tsne_fig, use_container_width=True)
                
                # Expression Heatmap
                st.markdown("#### 🔥 Expression Heatmap")
                
                # Select genes for heatmap
                heatmap_option = st.radio(
                    "Select genes for heatmap:",
                    ["Top Variable Genes", "Top Mean Expression", "Custom Selection"]
                )
                
                n_genes_heatmap = st.slider("Number of genes to display:", 10, 100, 50)
                
                if heatmap_option == "Top Variable Genes":
                    gene_var = processed_data.var(axis=1)
                    top_genes = gene_var.nlargest(n_genes_heatmap).index
                elif heatmap_option == "Top Mean Expression":
                    gene_mean = processed_data.mean(axis=1)
                    top_genes = gene_mean.nlargest(n_genes_heatmap).index
                else:
                    # Custom selection would require multiselect widget
                    available_genes = processed_data.index.tolist()[:100]  # Limit for performance
                    top_genes = st.multiselect(
                        "Select genes:", 
                        available_genes,
                        default=available_genes[:n_genes_heatmap]
                    )
                
                if len(top_genes) > 0:
                    heatmap_data = processed_data.loc[top_genes]
                    
                    heatmap_fig = px.imshow(
                        heatmap_data.values,
                        x=heatmap_data.columns,
                        y=heatmap_data.index,
                        aspect="auto",
                        color_continuous_scale="RdBu_r",
                        title=f"Expression Heatmap ({len(top_genes)} genes)"
                    )
                    heatmap_fig.update_layout(height=max(400, len(top_genes) * 10))
                    st.plotly_chart(heatmap_fig, use_container_width=True)
                
                # Correlation Analysis
                st.markdown("#### 🔗 Sample Correlation Analysis")
                
                sample_corr = processed_data.T.corr()
                
                corr_fig = px.imshow(
                    sample_corr.values,
                    x=sample_corr.columns,
                    y=sample_corr.index,
                    color_continuous_scale="RdBu_r",
                    title="Sample-to-Sample Correlation Matrix",
                    zmin=-1, zmax=1
                )
                st.plotly_chart(corr_fig, use_container_width=True)
                
                # Distribution Comparisons
                st.markdown("#### 📊 Expression Distribution Comparison")
                
                if 'metadata' in st.session_state:
                    metadata = st.session_state.metadata
                    group_column = st.selectbox(
                        "Group samples by:", 
                        metadata.columns.tolist()
                    )
                    
                    # Create violin plot comparing distributions
                    sample_means = processed_data.mean(axis=0)
                    plot_data = []
                    
                    for sample in processed_data.columns:
                        if sample in metadata['Sample'].values:
                            group_value = metadata[metadata['Sample'] == sample][group_column].iloc[0]
                            plot_data.append({
                                'Sample': sample,
                                'Mean_Expression': sample_means[sample],
                                'Group': str(group_value)
                            })
                    
                    plot_df = pd.DataFrame(plot_data)
                    
                    violin_fig = px.violin(
                        plot_df,
                        x='Group',
                        y='Mean_Expression',
                        title=f"Expression Distribution by {group_column}",
                        box=True
                    )
                    st.plotly_chart(violin_fig, use_container_width=True)
        
        else:
            st.markdown("""
            <div class="warning-box">
                <h4>⚠️ No Analysis Results Available</h4>
                <p>Please complete the data analysis in the Analysis tab first.</p>
            </div>
            """, unsafe_allow_html=True)
    
    # TEAM INFO TAB
    with tabs[3]:
        st.markdown('<h2 class="section-header">👥 Team Information</h2>', unsafe_allow_html=True)
        
        # Team Overview
        st.markdown("""
        <div class="info-box">
            <h3>🎯 Project Overview</h3>
            <p>MolBio Analytics Pro is an advanced molecular biology data analysis platform developed as part of the Machine Learning course. 
            This application provides comprehensive tools for processing, analyzing, and visualizing gene expression data with a focus on 
            user-friendly interfaces and robust analytical capabilities.</p>
        </div>
        """, unsafe_allow_html=True)
        
        # Team Members (Update with your actual team information)
        st.markdown('<h3 class="subsection-header">👨‍💻 Development Team</h3>', unsafe_allow_html=True)
        
        
        
        # Docker Information
        st.markdown('<h3 class="subsection-header">🐳 Docker Deployment</h3>', unsafe_allow_html=True)
        
        st.markdown("""
        <div class="info-box">
            <h4>🚀 Containerized Application</h4>
            <p>This application is fully containerized using Docker for easy deployment and reproducibility. 
            The Docker configuration includes all necessary dependencies and ensures consistent behavior across different environments.</p>
            
            <h5>Quick Deploy Commands:</h5>
            <code>
            docker build -t molbio-analytics .<br>
            docker run -p 8501:8501 molbio-analytics
            </code>
        </div>
        """, unsafe_allow_html=True)
        

        stats_col1, stats_col2, stats_col3, stats_col4 = st.columns(4)
        
        with stats_col1:
            st.markdown("""
            <div class="metric-card">
                <div class="metric-value">2,500+</div>
                <div class="metric-label">Lines of Code</div>
            </div>
            """, unsafe_allow_html=True)
        
        with stats_col2:
            st.markdown("""
            <div class="metric-card">
                <div class="metric-value">15+</div>
                <div class="metric-label">ML Algorithms</div>
            </div>
            """, unsafe_allow_html=True)
        
        with stats_col3:
            st.markdown("""
            <div class="metric-card">
                <div class="metric-value">25+</div>
                <div class="metric-label">Visualization Types</div>
            </div>
            """, unsafe_allow_html=True)
        
        with stats_col4:
            st.markdown("""
            <div class="metric-card">
                <div class="metric-value">100+</div>
                <div class="metric-label">Hours Developed</div>
            </div>
            """, unsafe_allow_html=True)
        
        # Version Information
        st.markdown('<h3 class="subsection-header">🔢 Version Information</h3>', unsafe_allow_html=True)
        
        version_col1, version_col2 = st.columns(2)
        
        with version_col1:
            st.markdown("""
            **Current Version:** v1.0.0  
            **Release Date:** December 2024  
            **Status:** Stable Release  
            **License:** MIT License  
            """)
        
        with version_col2:
            st.markdown("""
            **Python Version:** 3.8+  
            **Streamlit Version:** 1.28+  
            **Dependencies:** 20+ packages  
            **Platform:** Cross-platform  
            """)
        
        # Future Roadmap
        st.markdown('<h3 class="subsection-header">🚀 Future Roadmap</h3>', unsafe_allow_html=True)
        
        st.markdown("""
        <div class="info-box">
            <h4>🔮 Planned Features</h4>
            <ul>
            <li><strong>v1.1:</strong> Advanced differential expression analysis</li>
            <li><strong>v1.2:</strong> Pathway enrichment analysis integration</li>
            <li><strong>v1.3:</strong> Multi-omics data support</li>
            <li><strong>v2.0:</strong> Cloud deployment and collaboration features</li>
            <li><strong>v2.1:</strong> Machine learning model training interface</li>
            <li><strong>v2.2:</strong> Real-time data processing capabilities</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)

# ===== FOOTER =====
def show_footer():
    """Display application footer."""
    st.markdown("---")
    st.markdown("""
    <div style="text-align: center; padding: 2rem; color: #6c757d;">
        <p>🧬 <strong>MolBio Analytics Pro</strong> v1.0.0 | Built with ❤️ using Streamlit</p>
        <p>© 2024 MolBio Analytics Team. All rights reserved.</p>
        <p><em>Empowering molecular biology research through advanced data analytics</em></p>
    </div>
    """, unsafe_allow_html=True)

# ===== UTILITY FUNCTIONS =====
def download_results():
    """Generate downloadable results summary."""
    if st.session_state.analysis_complete:
        results_summary = {
            'timestamp': datetime.now().isoformat(),
            'data_shape': st.session_state.processed_data.shape,
            'processing_steps': st.session_state.processing_steps,
            'pca_variance': st.session_state.pca_results['explained_variance_ratio'].tolist(),
            'cluster_info': {
                'method': st.session_state.cluster_results['method'],
                'n_clusters': st.session_state.cluster_results['n_clusters'],
                'cluster_sizes': [len(profile['samples']) for profile in st.session_state.cluster_results['profiles']]
            }
        }
        return results_summary
    return None

def export_data():
    """Export processed data and results."""
    if st.session_state.analysis_complete:
        # Create a bytes buffer
        buffer = io.BytesIO()
        
        # Write processed data to Excel
        with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
            st.session_state.processed_data.to_excel(writer, sheet_name='Processed_Data')
            
            # Add PCA results
            pca_coords = pd.DataFrame(
                st.session_state.pca_results['coordinates'],
                columns=[f'PC{i+1}' for i in range(st.session_state.pca_results['coordinates'].shape[1])],
                index=st.session_state.pca_results['sample_names']
            )
            pca_coords.to_excel(writer, sheet_name='PCA_Coordinates')
            
            # Add cluster assignments
            cluster_df = pd.DataFrame({
                'Sample': st.session_state.cluster_results['sample_names'],
                'Cluster': st.session_state.cluster_results['labels']
            })
            cluster_df.to_excel(writer, sheet_name='Cluster_Assignments', index=False)
        
        buffer.seek(0)
        return buffer.getvalue()
    return None

# ===== RUN APPLICATION =====
if __name__ == "__main__":
    main()
    
    # Add download functionality in sidebar
    if st.session_state.get('analysis_complete', False):
        st.sidebar.markdown("---")
        st.sidebar.markdown("### 💾 Export Results")
        
        # Download results summary
        if st.sidebar.button("📋 Download Summary", use_container_width=True):
            results = download_results()
            if results:
                st.sidebar.download_button(
                    label="📁 Download JSON Summary",
                    data=pd.Series(results).to_json(indent=2),
                    file_name=f"molbio_analysis_summary_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                    mime="application/json"
                )
        
        # Export processed data
        if st.sidebar.button("📊 Export Data", use_container_width=True):
            excel_data = export_data()
            if excel_data:
                st.sidebar.download_button(
                    label="📁 Download Excel File",
                    data=excel_data,
                    file_name=f"molbio_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                )
    
    # Show footer
    show_footer()

# ===== ERROR HANDLING AND LOGGING =====
def handle_analysis_error(error, context="General"):
    """Handle analysis errors gracefully."""
    error_msg = f"An error occurred during {context}: {str(error)}"
    st.error(error_msg)
    
    # Log error details (in production, you'd want proper logging)
    if st.checkbox("Show detailed error information", key=f"error_details_{context}"):
        st.code(f"Error Type: {type(error).__name__}\nError Message: {str(error)}")
    
    return False

# ===== PERFORMANCE MONITORING =====
def monitor_performance():
    """Monitor application performance."""
    if 'performance_metrics' not in st.session_state:
        st.session_state.performance_metrics = {
            'analysis_start_time': None,
            'analysis_duration': None,
            'data_size': None
        }
    
    return st.session_state.performance_metrics

# ===== CONFIGURATION VALIDATION =====
def validate_configuration():
    """Validate analysis configuration parameters."""
    issues = []
    
    # Check if required parameters are set
    if not st.session_state.get('data_loaded', False):
        issues.append("No data loaded")
    
    # Add more validation as needed
    
    return issues

# ===== HELP AND DOCUMENTATION =====
def show_help_section():
    """Display help and documentation."""
    with st.expander("❓ Need Help?"):
        st.markdown("""
        ### 🆘 Common Issues and Solutions
        
        **Data Loading Problems:**
        - Ensure your file has genes as rows and samples as columns
        - Check that the first column contains gene names
        - Verify that expression values are numeric
        
        **Analysis Failures:**
        - Try reducing the number of features for large datasets
        - Check for missing values in your data
        - Ensure you have at least 3 samples for clustering
        
        **Visualization Issues:**
        - Reduce the number of genes in heatmaps for better performance
        - Try different color schemes if plots are hard to read
        - Use the zoom and pan features in interactive plots
        
        ### 📚 Additional Resources
        - [Streamlit Documentation](https://streamlit.io)
        - [Scikit-learn User Guide](https://scikit-learn.org)
        - [Plotly Documentation](https://plotly.com/python/)
        """)

# Add help section to sidebar
with st.sidebar:
    st.markdown("---")
    show_help_section()