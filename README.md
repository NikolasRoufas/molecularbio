# 🧬 MolBio Analytics Pro

<div align="center">

![Python](https://img.shields.io/badge/python-v3.8+-blue.svg)
![Streamlit](https://img.shields.io/badge/streamlit-v1.28+-red.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
![Status](https://img.shields.io/badge/status-stable-brightgreen.svg)

**Advanced Molecular Biology Data Analysis Platform**

*Empowering molecular biology research through advanced data analytics*
**created for Ionian University**

</div>

## 📋 Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage Guide](#usage-guide)
- [Docker Deployment](#docker-deployment)
- [API Reference](#api-reference)
- [Data Formats](#data-formats)
- [Analysis Methods](#analysis-methods)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)

## 🎯 Overview

MolBio Analytics Pro is a comprehensive web-based platform designed for analyzing molecular biology data, with a primary focus on gene expression analysis. Built with Streamlit and powered by advanced machine learning algorithms, it provides researchers with an intuitive interface for processing, analyzing, and visualizing complex biological datasets.

### Key Capabilities

- **Data Processing**: Advanced preprocessing pipelines with filtering, normalization, and transformation
- **Quality Control**: Comprehensive QC metrics and visualizations
- **Dimensionality Reduction**: PCA and t-SNE analysis with interactive visualizations
- **Clustering Analysis**: K-means and hierarchical clustering with optimization
- **Interactive Visualizations**: Rich, interactive plots using Plotly
- **Export Functionality**: Download results in multiple formats

## ✨ Features

### 🔧 Data Processing
- **Multi-format Support**: CSV, Excel (XLSX), and TSV files
- **Smart Data Cleaning**: Automatic handling of missing values and data type conversion
- **Advanced Filtering**: Remove low-expression genes and high-zero genes
- **Multiple Transformations**: Log2, Log10, and natural log transformations
- **Normalization Methods**: TPM and quantile normalization
- **Scaling Options**: Standard and robust scaling

### 📊 Quality Control
- **Comprehensive Metrics**: Library sizes, detection rates, expression distributions
- **Visual QC Dashboard**: Interactive plots for data quality assessment
- **Sample-wise Statistics**: Individual sample quality metrics
- **Gene-wise Statistics**: Per-gene expression characteristics

### 🔬 Advanced Analytics
- **Principal Component Analysis**: Variance explained, loadings analysis, feature importance
- **Clustering Analysis**: K-means with elbow method optimization, hierarchical clustering
- **t-SNE Visualization**: Non-linear dimensionality reduction
- **Statistical Analysis**: Distribution analysis, correlation matrices

### 📈 Visualizations
- **Interactive Plots**: Zoom, pan, and hover functionality
- **Expression Heatmaps**: Customizable gene selection and clustering
- **PCA Plots**: Color by metadata, variance explained plots
- **Cluster Visualizations**: 2D projections with cluster assignments
- **Correlation Matrices**: Sample-to-sample relationships

## 🚀 Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager

### Local Installation

1. **Clone the repository**


2. **Create virtual environment**
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. **Install dependencies**
```bash
pip install -r requirements.txt
```

4. **Run the application**
```bash
streamlit run app.py
```

### Dependencies

```
streamlit>=1.28.0
pandas>=1.5.0
numpy>=1.21.0
plotly>=5.0.0
scikit-learn>=1.1.0
seaborn>=0.11.0
scipy>=1.7.0
openpyxl>=3.0.0
```

## 🎯 Quick Start

### 1. Launch the Application
```bash
streamlit run app.py
```

### 2. Load Your Data
- **Option A**: Upload your gene expression file (CSV, XLSX, or TSV)
- **Option B**: Generate a sample dataset for testing

### 3. Configure Analysis Parameters
- Set preprocessing options in the sidebar
- Choose analysis methods (PCA, clustering, t-SNE)
- Adjust visualization parameters

### 4. Run Analysis
- Navigate to the "Analysis" tab
- Click "Process Data" to start the analysis
- Monitor progress and view processing summary

### 5. Explore Results
- Navigate to the "Results" tab
- Explore quality control metrics
- Visualize PCA and clustering results
- Generate advanced plots and export data

## 📖 Usage Guide

### Data Input Requirements

Your input data should be structured as follows:
- **Rows**: Genes/Features
- **Columns**: Samples
- **Values**: Expression levels (numeric)
- **Index**: Gene names/IDs
- **Headers**: Sample names/IDs

Example structure:
```
Gene_ID    Sample1  Sample2  Sample3  ...
GENE001    12.45    8.32     15.67    ...
GENE002    3.21     9.87     4.56     ...
GENE003    45.12    23.45    38.90    ...
...
```

### Analysis Workflow

#### Step 1: Data Preprocessing
1. **Filter Low-Quality Data**
   - Remove genes with high zero fractions
   - Filter low-expression genes
   
2. **Transform Data**
   - Apply log transformation (recommended)
   - Choose appropriate normalization method
   
3. **Scale Data** (optional)
   - Standard scaling for PCA
   - Robust scaling for outlier-resistant analysis

#### Step 2: Quality Control
- Review basic statistics
- Examine expression distributions
- Check sample correlations
- Identify potential outliers

#### Step 3: Dimensionality Reduction
- **PCA Analysis**
  - Determine optimal number of components
  - Analyze variance explained
  - Identify top contributing genes
  
- **t-SNE Analysis** (optional)
  - Non-linear dimensionality reduction
  - Adjust perplexity and iterations

#### Step 4: Clustering Analysis
- **K-means Clustering**
  - Automatic optimal K determination
  - Silhouette score optimization
  
- **Hierarchical Clustering**
  - Ward linkage method
  - Dendrogram visualization

#### Step 5: Results Interpretation
- Analyze cluster characteristics
- Identify marker genes
- Export results for further analysis

### Advanced Configuration

#### Preprocessing Parameters
```python
processing_params = {
    'filter_zeros': True,
    'zero_threshold': 0.5,
    'filter_low_expr': True,
    'min_expression': 1.0,
    'log_transform': True,
    'log_type': 'log2',
    'normalize': False,
    'norm_method': 'tpm',
    'scale_data': False,
    'scale_method': 'standard'
}
```

#### Analysis Parameters
```python
analysis_params = {
    'n_pca_components': 10,
    'clustering_method': 'kmeans',
    'n_clusters': 0,  # 0 for automatic
    'run_tsne': False,
    'tsne_perplexity': 30,
    'tsne_iterations': 1000
}
```

## 🐳 Docker Deployment

### Build Docker Image
```bash
docker build -t molbio-analytics-pro .
```

### Run Container
```bash
docker run -p 8501:8501 molbio-analytics-pro
```

### Docker Compose (Recommended)
```yaml
version: '3.8'
services:
  molbio-analytics:
    build: .
    ports:
      - "8501:8501"
    volumes:
      - ./data:/app/data
    environment:
      - STREAMLIT_SERVER_PORT=8501
      - STREAMLIT_SERVER_ADDRESS=0.0.0.0
```

### Dockerfile
```dockerfile
FROM python:3.9-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

EXPOSE 8501

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

ENTRYPOINT ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]
```

## 🔧 API Reference

### Core Classes

#### `DataProcessor`
```python
@staticmethod
def generate_realistic_sample_data() -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
    # Generates biologically realistic sample data

@staticmethod
def preprocess_data(data: pd.DataFrame, **kwargs) -> Tuple[pd.DataFrame, List[str]]
    # Preprocesses molecular biology data
```

#### `QualityControl`
```python
@staticmethod
def generate_qc_metrics(data: pd.DataFrame, metadata: pd.DataFrame = None) -> Dict
    # Generates comprehensive QC metrics
```

#### `AdvancedAnalysis`
```python
@staticmethod
def perform_pca(data: pd.DataFrame, n_components: int = 10, metadata: pd.DataFrame = None) -> Dict
    # Performs PCA analysis

@staticmethod
def perform_clustering(data: pd.DataFrame, method: str = 'kmeans', n_clusters: int = None, **kwargs) -> Dict
    # Performs clustering analysis
```

#### `VisualizationEngine`
```python
@staticmethod
def create_qc_dashboard(data: pd.DataFrame, qc_metrics: Dict) -> go.Figure
    # Creates QC dashboard

@staticmethod
def create_pca_visualization(pca_results: Dict, metadata: pd.DataFrame = None, color_by: str = None) -> go.Figure
    # Creates PCA visualization

@staticmethod
def create_clustering_visualization(cluster_results: Dict, pca_results: Dict = None) -> go.Figure
    # Creates clustering visualization
```

## 📁 Data Formats

### Supported Input Formats

1. **CSV Files (.csv)**
   - Comma-separated values
   - UTF-8 encoding recommended
   - First column as gene index

2. **Excel Files (.xlsx)**
   - Microsoft Excel format
   - Single worksheet support
   - First column as gene index

3. **TSV Files (.tsv)**
   - Tab-separated values
   - UTF-8 encoding recommended
   - First column as gene index

### Output Formats

1. **JSON Summary**
   - Analysis parameters
   - Key results and metrics
   - Timestamp and metadata

2. **Excel Results**
   - Processed data
   - PCA coordinates
   - Cluster assignments
   - Multiple worksheets

3. **Interactive HTML**
   - Plotly visualizations
   - Embeddable in reports
   - Full interactivity preserved

## 🧪 Analysis Methods

### Principal Component Analysis (PCA)
- **Purpose**: Dimensionality reduction and feature extraction
- **Method**: Singular Value Decomposition
- **Scaling**: Automatic standardization
- **Output**: Component coordinates, explained variance, loadings

### Clustering Analysis

#### K-means Clustering
- **Algorithm**: Lloyd's algorithm with k-means++
- **Optimization**: Elbow method and silhouette analysis
- **Initialization**: 10 random initializations
- **Distance Metric**: Euclidean distance

#### Hierarchical Clustering
- **Linkage**: Ward method (minimizes within-cluster variance)
- **Distance Metric**: Euclidean distance
- **Output**: Dendrogram and cluster assignments

### t-SNE (t-Distributed Stochastic Neighbor Embedding)
- **Purpose**: Non-linear dimensionality reduction
- **Perplexity**: Adjustable (5-50 range)
- **Iterations**: Configurable (250-2000)
- **Initialization**: Random

### Quality Control Metrics

#### Sample-Level Metrics
- Total expression counts
- Number of detected genes
- Mean expression level
- Coefficient of variation

#### Gene-Level Metrics
- Mean expression across samples
- Expression variance
- Detection rate (fraction of samples with expression > 0)
- Maximum expression level

## 🐛 Troubleshooting

### Common Issues

#### Data Loading Problems
**Problem**: File upload fails
```
Solution:
- Check file format (CSV, XLSX, TSV only)
- Ensure first column contains gene names
- Verify numeric expression values
- Check file encoding (UTF-8 recommended)
```

**Problem**: Memory error with large files
```
Solution:
- Filter data before upload
- Use sample dataset for testing
- Increase system memory
- Process data in chunks
```

#### Analysis Failures
**Problem**: PCA analysis fails
```
Solution:
- Check for missing values
- Ensure minimum 2 samples and 2 genes
- Verify numeric data types
- Try reducing number of components
```

**Problem**: Clustering optimization fails
```
Solution:
- Ensure minimum 3 samples
- Check for infinite/NaN values
- Try manual cluster number selection
- Use hierarchical clustering instead
```

#### Visualization Issues
**Problem**: Plots don't render
```
Solution:
- Check browser compatibility
- Clear browser cache
- Disable ad blockers
- Try different browser
```

**Problem**: Heatmap too large/slow
```
Solution:
- Reduce number of genes displayed
- Filter top variable genes
- Use gene selection options
- Increase browser memory
```

### Error Messages

| Error | Cause | Solution |
|-------|-------|----------|
| `TypeError: perform_pca() got unexpected keyword` | Incorrect parameter passing | Remove unsupported parameters |
| `ValueError: Input contains NaN` | Missing values in data | Enable data cleaning options |
| `MemoryError` | Dataset too large | Filter data or use sampling |
| `KeyError: 'column_name'` | Missing metadata column | Check metadata file structure |

### Performance Optimization

#### For Large Datasets
1. **Filter aggressively**: Remove low-expression genes
2. **Sample data**: Use representative subset
3. **Reduce features**: Focus on most variable genes
4. **Optimize parameters**: Lower PCA components, fewer clusters

#### Memory Management
1. **Close unused tabs**: Free browser memory
2. **Restart application**: Clear session state
3. **Use chunking**: Process data in smaller batches
4. **Monitor resources**: Check system memory usage

## 🤝 Contributing

We welcome contributions to MolBio Analytics Pro! Here's how you can help:

### Development Setup
1. Fork the repository
2. Create a feature branch
3. Install development dependencies
4. Make your changes
5. Run tests
6. Submit a pull request

### Contribution Guidelines
- Follow PEP 8 style guidelines
- Add docstrings to new functions
- Include unit tests for new features
- Update documentation as needed
- Test with sample datasets

### Areas for Contribution
- **New Analysis Methods**: Implement additional statistical tests
- **Visualization Enhancements**: Add new plot types
- **Performance Optimization**: Improve processing speed
- **Documentation**: Expand user guides and tutorials
- **Testing**: Add comprehensive test coverage

## 📊 Project Statistics

- **Lines of Code**: 2,500+
- **Functions**: 50+
- **Classes**: 4 main classes
- **Dependencies**: 20+ packages
- **Visualization Types**: 25+
- **Analysis Methods**: 15+

## 🗺️ Roadmap

### Version 1.1 (Next Release)
- [ ] Differential expression analysis
- [ ] Gene set enrichment analysis
- [ ] Batch effect correction
- [ ] Enhanced metadata support

### Version 1.2
- [ ] Pathway analysis integration
- [ ] Multi-omics data support
- [ ] Advanced statistical tests
- [ ] Custom gene signatures

### Version 2.0 (Future)
- [ ] Cloud deployment options
- [ ] Collaboration features
- [ ] API endpoints
- [ ] Real-time processing
- [ ] Machine learning model training

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

```
MIT License

Copyright (c) 2025 Nikolaos Roufas

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```



## 🏆 Acknowledgments

- Built with [Streamlit](https://streamlit.io/)
- Powered by [scikit-learn](https://scikit-learn.org/)
- Visualizations by [Plotly](https://plotly.com/)
- Data processing with [pandas](https://pandas.pydata.org/)
- Statistical analysis with [SciPy](https://scipy.org/)

---

<div align="center">


</div>
