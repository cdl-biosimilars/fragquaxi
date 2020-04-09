#' Amino acid residues.
#'
#' A dataset containing abbreviations and formulas
#' for the 22 proteinogenic amino acids.
#'
#' @format A data frame with 22 rows and 4 variables:
#' \describe{
#' \item{abbreviation}{single-letter abbreviation}
#' \item{three_letter}{three-letter abbreviation}
#' \item{name}{name}
#' \item{formula}{formula for the residue (i.e., for an amino acid
#'                incorporated into a polypeptide chain)}
#' }
"amino_acids"

#' Average and monoisotopic atomic masses.
#'
#' A dataset containing average and monoisotopic masses
#' for the chemical elements. Masses are missing in case of elements
#' for which no stable isotope is known.
#'
#' @format A data frame with 118 rows and 6 variables:
#' \describe{
#' \item{element}{chemical symbol}
#' \item{name}{IUPAC name}
#' \item{Z}{atomic number (number of protons in the nucleus)}
#' \item{average}{average atomic mass as defined by IUPAC}
#' \item{A}{mass number (total number of protons and neutrons in the nucleus)
#'          of the most abundant isotope}
#' \item{monoisotopic}{atomic mass of the most abundant isotope}
#' }
#'
#' @source
#' * Average masses: Atomic weights of the elements 2013, Tables 2 and 3
#'   (<https://doi.org/10.1515/pac-2015-0305>)
#' * Isotope abundances: Isotopic compositions of the elements 2013, Table 1
#'   (<https://doi.org/10.1515/pac-2015-0503>)
#' * Monoisotopic masses: The AME16 Atomic Mass Evaluation
#'   (<https://www-nds.iaea.org/amdc/ame2016/mass16.txt>)
"atomic_masses"

#' Common monosaccharide residues.
#'
#' A dataset describing molecular compositions for common monosaccharide
#' residues.
#'
#' @format A data frame with 83 rows and 5 variables:
#' \describe{
#' \item{abbreviation}{abbreviation}
#' \item{type}{carbohydrate type}
#' \item{formula}{formula for the residue (i.e., for a monosaccharide
#'                incorporated into a glycan)}
#' \item{short_name}{short name}
#' \item{systematic_name}{IUPAC name}
#' }
#'
#' @source
#' * Abbreviations and names: Symbol Nomenclature for Glycans (SNFG),
#'   2019 edition, Table 3 (<https://www.ncbi.nlm.nih.gov/glycans/snfg.html>)
#' * Formulas: PubChem (<https://pubchem.ncbi.nlm.nih.gov/>)
"monosaccharides"

#' Common post-translational modifications.
#'
#' A dataset describing molecular compositions for common post-translational
#' modifications.
#'
#' @format A data frame with 72 rows and 3 variables:
#' \describe{
#' \item{abbreviation}{abbreviation}
#' \item{name}{description}
#' \item{formula}{correction formula}
#' }
#'
#' @source
#' * Abbreviations and name: ExPASy FindMod tool
#'   (<https://web.expasy.org/findmod/findmod_masses.html>)
#' * Correction formulas: UniProt Controlled vocabulary of
#'   posttranslational modifications (PTM), release 2020_01 of 26-Feb-2020
#'   (<https://www.uniprot.org/docs/ptmlist>)
"ptms"

#' Exemplary PTM compositions ("modcoms").
#'
#' A dataset describing PTM compositions.
#'
#' @format A data frame with 6 rows and 4 variables:
#' \describe{
#' \item{modcom_name}{abbreviation}
#' \item{Hex}{number of hexoses}
#' \item{HexNAc}{number of N-acetylhexosamines}
#' \item{Fuc}{number of fucoses}
#' }
#'
"sample_modcoms"
