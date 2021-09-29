\name{wheat}
\alias{wheatHTP}
\alias{Y}
\alias{X}
\alias{M}
\alias{CV}
\alias{genCOV_xy}
\alias{genCOV_yy}
\alias{resCOV_yy}
\alias{wheatHTP.E1}
\alias{wheatHTP.E2}
\alias{wheatHTP.E3}
\alias{wheatHTP.E4}
\title{Wheat dataset}
\description{

The dataset consists of 1,092 inbred wheat lines grouped into 39 trials and grown during the 2013-2014 season at the Norman Borlaug experimental research station in Ciudad Obregon, Sonora, Mexico. Each trial consisted of 28 breeding lines that were arranged in an alpha-lattice design with three replicates and six sub-blocks. The trials were grown in four different environments:
\itemize{
\item E1: Flat-Drought (sowing in flat with irrigation of 180 mm through drip system)
\item E2: Bed-2IR (sowing in bed with 2 irrigations approximately 250 mm)
\item E3: Bed-5IR (bed sowing with 5 normal irrigations)
\item E4: Bed-EHeat (bed sowing 30 days before optimal planting date with 5 normal irrigations approximately 500 mm)
}

\ifelse{html}{\out{<b>1. Phenotypic data.</b>}}{\eqn{\textbf{1. Phenotypic data.}}{1. Phenotypic data.}}
Measurements of grain yield (YLD) were reported as the total plot yield after maturity. Records for YLD are reported as adjusted means from which trial, replicate and sub-block effects were removed. Measurements for days to heading (DTH), days to maturity (DTM), and plant height (PH) were recorded only in the first replicate at each trial and thus no phenotype adjustment was made.

\ifelse{html}{\out{<b>2. Reflectance data.</b>}}{\eqn{\textbf{2. Reflectance data.}}{2. Reflectance data.}}
Reflectance data was collected from the fields using both infrared and hyper-spectral cameras mounted on an aircraft on 9 different dates (time-points) between January 10 and March 27th, 2014. During each flight, data from 250 wavelengths ranging from 392 to 850 nm were collected for each pixel in the pictures. The average reflectance of all the pixels for each wavelength was calculated from each of the geo-referenced trial plots and reported as each line reflectance. Data for reflectance and Green NDVI and Red NDVI are reported as adjusted phenotypes from which trial, replicate and sub-block effects were removed. Each data-point matches to each data-point in phenotypic data.

\ifelse{html}{\out{<b>3. Marker data.</b>}}{\eqn{\textbf{3. Marker data.}}{3. Marker data.}} 
Lines were sequenced for GBS at 192-plexing on Illumina HiSeq2000 or HiSeq2500 with 1 x 100 bp reads. SNPs were called across all lines anchored to the genome assembly of Chinese Spring (International Wheat Genome Sequencing Consortium 2014). Next, SNP were extracted and filtered so that lines >50\% missing data were removed. Markers were recoded as –1, 0, and 1, corresponding to homozygous for the minor allele, heterozygous, and homozygous for the major allele, respectively. Next, markers with a minor allele frequency <0.05 and >15\% of missing data were removed. Remaining SNPs with missing values were imputed using the mean of the observed marker genotypes at a given locus.

\ifelse{html}{\out{<b>Replicated and Un-replicated data.</b>}}{\eqn{\textbf{Replicated and Un-replicated data.}}{Replicated and Un-replicated data.}}
The CRAN version includes the \code{wheatHTP} dataset containing (un-replicated) phenotypic and reflectance (latest time-point only) data from the environment E3 (Bed-5IR) only. Marker data is also included in the dataset. The phenotypic and reflectance data are averages (from adjusted phenotypes) for 967 lines evaluated in 35 trials for which marker information on 3,274 SNPs is available.

The GitHub (development) version of the SFSI R-package (https://github.com/MarcooLopez/SFSI) includes also the \code{wheatHTP} dataset plus \code{wheatHTP.E1},...,\code{wheatHTP.E4} datasets containing replicated (adjusted) phenotypic and reflectance data for all four environments, respectively.

\ifelse{html}{\out{<b>Cross-validation partitions.</b>}}{\eqn{\textbf{Cross-validation partitions.}}{Cross-validation partitions.}}
Random partitions of training-testing sets are provided for the un-replicated data (\code{wheatHTP} dataset) containing 967 individuals and 35 trials. Data from 10 entire trials (30\% of 35 the trials) were randomly assigned to the testing set, and data from the remaining 25 trials were assigned to the training set. Partitions consist of arrays of 967 indices with 1's (training) and 2's (testing). 10 different partitions were generated and storaged in a matrix of 967 rows and 10 columns. 

\ifelse{html}{\out{<b>Genetic covariances.</b>}}{\eqn{\textbf{Genetic covariances.}}{Genetic covariances.}}
Multi-variate Gaussian mixed models were fitted to phenotypes from the un-replicated data (\code{wheatHTP} dataset) containing 967 individuals. Bi-variate models were fitted to YLD with each of the 250 wavelengths. Tetra-variate models were fitted for YLD, DTH, DTM, and PH traits. All models were fitted within training set to scaled (null mean and unit variance) phenotypes for all partitions (provided matrix). Bayesian models were implemented using the 'Multitrait' function from the \code{BGLR} R-package with 25,000 iterations discarding 5,000 runs for burning. A marker-derived relationships matrix as in VanRaden (2008) was used to model between-individuals genetic effects; both between-traits genetic and residual covariances were assumed unstructured. 

Genetic covariances between YLD and each wavelength are storaged in a matrix of 250 rows and 10 columns (partitions). Genetic and residual covariances matrices among YLD, DTH, DTM, and PH are storaged in a list with 10 elements (partitions). 
}

\usage{
  data(wheatHTP)
  # data(wheatHTP.E1) # GitHub version
  # data(wheatHTP.E2) # GitHub version
  # data(wheatHTP.E3) # GitHub version
  # data(wheatHTP.E4) # GitHub version
}

\format{
Replicated data (\code{wheatHTP.E1},...,\code{wheatHTP.E4} datasets):
\itemize{
\item \code{Y}: (matrix) phenotypic data for YLD, DTH, DTM, and PH; and the trial in which each genotype was tested. 
\item \code{WL}: (9-dimensional list) reflectance data for time-points 1,2,...,9.
\item \code{VI}: (9-dimensional list) green and red NDVI for time-points 1,2,...,9. 
}

Un-replicted data (\code{wheatHTP} dataset):
\itemize{
\item \code{Y}: (matrix) phenotypic data for YLD, DTH, DTM, and PH; and the trial in which each genotype was tested. 
\item \code{M}: (matrix) marker data with SNPs in columns.
\item \code{X}: (matrix) reflectance data for time-point 9.
\item \code{CV}: (matrix) array of 1's (training) and  2's (testing) for 10 partitions (in columns). 
\item \code{genCOV_xy}: (matrix) genetic covariances between YLD and each reflectance trait for each partition (in columns). 
\item \code{genCOV_yy}: (10-dimensional list) genetic covariances matrices among YLD, DTH, DTM, and PH for partitions 1,2,...,10.
\item \code{resCOV_yy}: (10-dimensional list) residual covariances matrices among YLD, DTH, DTM, and PH for partitions 1,2,...,10.
}
}
\references{
Perez-Rodriguez P, de los Campos G (2014). Genome-wide regression and prediction with the BGLR statistical package. \emph{Genetics}, \bold{198}, 483–495.

VanRaden PM (2008). Efficient methods to compute genomic predictions. \emph{Journal of Dairy Science}, \bold{91}(11), 4414–4423.
}
\source{
  International Maize and Wheat Improvement Center (CIMMYT), Mexico.
}

\keyword{datasets}