/**
 * This package contains classes to set up and run an NPAIRS analysis; look here for an 
 * overview of NPAIRS analysis technique.
 *
 * <p></p>
 * <h1><a name="Getting Started"></a>Getting Started</h1>
 * <ol>
 *  <li><a href="#Introduction">Introduction</a> </li> 
 *  <li><a href="#Overview">Overview</a></li>
 *  <li><a href="#Details">Details</a></li>
 * </ol>
 * <h2><a name="Introduction"></a>1. Introduction</h2>
 * 
 * <blockquote>NPAIRS provides a resampling framework for combining prediction metrics
with the reproducibility of the brain-activation patterns, or statistical parametric
maps (SPM), as a data-driven substitute for ROCs. However, any
measure of similarity between patterns extracted from independent data sets
is subject to an unknown bias (Afshinpour et al., in press). To obtain combined
prediction and reproducibility values Strother et al. (2002; Kjems et al.,
2002) proposed a novel split-half resampling framework dubbed NPAIRS and
applied it first to PET and later to fMRI (see Strother et al., 2004; LaConte
et al., 2003; Yourganov et al., in press). While NPAIRS may be applied to
any analysis model we have focused on LDA built on a regularized PCA basis
(i.e., PDA). This allows us to (1) regularize the model by choosing soft (e.g.,
ridge) or hard thresholds on the PCA eigenspectrum or other basis set (e.g.,
tensor product splines) (2) maintain the link to covariance decomposition
previously used with PET for elucidating network structures, and (3) easily
produce robust whole-brain activation maps useful for discovering features of
brain function and/or disease. (Strother et al., <i>Comp. Stat. in 
Neuroimaging</i>, 2010)</blockquote>
 * </p> 
 * <p> Here is a diagram outlining the NPAIRS method :</p>
 * <p><img src="doc-files/NPAIRS_diagram_Strother2004.png" 
 * alt="NPAIRS diagram" width="70%">
 * (Image (c) Strother, S., 2004)
 * 
 * <h2><a name="Overview"></a>2. Overview</h2>
 * Here is a diagram showing an overview of the workflow in an NPAIRS analysis:
 * <p><img src="doc-files/NPAIRS_overview_workflow_big.png"
 * width="100%"></p>
 * <p><b><sup>1</sup><i>#PCs fed into CVA<sup>2</sup> = p<sub>i</sub> in ith analysis</i></b>
 * <p><b><sup>2</sup><i>See <a href="#PCA+CVA">PCA+CVA</a> flowchart</i></b>
 * <p><b><sup>3</sup><i>See <a href="#Splits analysis">Splits analysis</a> flowchart</i></b> 
 * <p><b><sup>4</sup><i>See <a href="#Procrustes">Procrustes</a> reference-matching flowchart</i></b>
 * 
 * <p><a href="#Introduction">Back</a> to Introduction</p>
 * 
 * <h2><a name="Details"></a>3. Details</h2>
 * <p>Here are details of the steps in an NPAIRS analysis.
 * <p><a name="PCA+CVA"></a><b>1. PCA + CVA</b></h2>
 * <p><img src="doc-files/PCA_CVA_big.png"
 * width="100%"></p>
 * 
 * <p><a href="#Overview">Back</a> to Overview</p>
 * 
 * <p><a name="Splits analysis"</a><b>2. Analysis of data splits</b></h2>
 * <p><img src="doc-files/Splits_analysis_big.png"
 * width="100%"></p> 
 * 
 * <p><b><sup>1</sup><i>See <a href="#PCA+CVA">PCA+CVA</a> flowchart</i></b>
 * <p><b><sup>2</sup><i>See <a href="#Procrustes">Procrustes</a> reference-matching flowchart</i></b>
 *
 * <p><a href="#Overview">Back</a> to Overview</p>
 * 
 * <p><a name="Procrustes"</a><b>3. Matching splits results to reference full-data results with Procrustes</b></h2>
 * <p><img src="doc-files/Procrust_big.png"
 * width="100%"></p>
 * 
 * <p><a href="#Overview">Back</a> to Overview
 * 
 *
 */
package npairs;