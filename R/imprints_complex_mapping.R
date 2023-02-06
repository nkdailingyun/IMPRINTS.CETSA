#' imprints_complex_mapping
#'
#' Function to map proteins to known protein complex database
#'
#' @param data dataset after calculating the relative protein abundance differences
#' and deriving the average of the profile, i.e., imprints_average(), make sure the columns
#' with readings are named in the format like "37C_Treatment"
#' @param categorytable dataset after imprints_score() or imprints_globalview(), the data
#' should contain one column named after "category", to allow focusing on a subset of the data,
#' say the "Changed" proteins
#' @param set a single character to indicate the sample name to analyze if any
#' @param treatment a single character to indicate the sample name to analyze
#' @param targetcategory the default value is c("NC","CN","CC")
#' @param complexdatabase by default is core Corum database, 03.09.2018 Corum 3.0 current release
#' @param organism choose from Human/Mouse/Rat, default is Human
#' @param complexID a vector of complexID to retrieve from the data, if specified
#' @param minsubunitsIdentified the minimal number of subunits should be identified,
#' default value is 3
#' @param cutisoformbycount when there are several isoforms mapped in a complex, by default to
#' filter away the isoforms with comparatively fewer count number
#'
#'
#' @import dplyr
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'   MOLM <- imprints_complex_mapping(MOLM, MOLM_score, treatment="TNF", targetcategory=c("NC","CC"))
#'   MOLM <- imprints_complex_mapping(MOLM, MOLM_category, treatment="TNF", complexID=c(1,10,1234))
#' }
#'
#'

imprints_complex_mapping <- function(data, categorytable=NULL, set=NULL, treatment=NULL,
                                     targetcategory=c("NC","CN","CC"),
                                     complexdatabase="Corum", organism="Human",complexID=NULL,
                                     minsubunitsIdentified=3, cutisoformbycount=TRUE) {

  if (complexdatabase!="Corum") { stop("Only core Corum database is supported for now") }
  if (!organism %in% c("Human","Mouse","Rat")) { stop("Provide the right organism") }
  complexdb <- subset(IMPRINTS.CETSA::coreCorum, Organism==organism)
  if (length(complexID)) {
    if (class(complexID)=="character") { complexID <- as.numeric(complexID) }
    complexdb <- subset(complexdb, ComplexID %in% complexID)
    if (nrow(complexdb)==0) { stop("Make sure the provided complexID is present in core Corum database") }
  }

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  if (length(treatment)==1) {
    treat=treatment
    print(paste0("To check the data under the treatment of ", treat))
    data <- subset(data, treatment==treat)
    sel <- grep(treatment,names(data))
    if (length(set)==1) {
      print(paste0("To check the data within the set of ", set))
      sel1 <- grep(set,names(data))
      sel <- intersect(sel,sel1)
    }
    if (length(sel)) {
      data <- data[ ,c(1,2,sel,(ncol(data)-2):ncol(data))]
      # print(head(data))
    }
  } else {
    stop("Need to specify a single character to indicate the condition name")
  }

  if (length(categorytable)>0 & length(grep("category",names(categorytable)))) {
    if ( length(targetcategory)>0 ) {
      categorytable <- subset(categorytable, category %in% targetcategory)
    }
    categorytable <- subset(categorytable, treatment==treat)
    if (length(set) & length(grep("category",names(categorytable)))) {
      setname = set
      categorytable <- subset(categorytable, set==setname)
    }
    data <- merge(data, categorytable[,c("id","category")], by="id")
  }

  if (nrow(data)!=0) {
    print(paste0(nrow(data), " protein entries to be assigned to Corum database..."))
  } else {
    stop("Make sure the right categories were specified.")
  }
  # print(head(data))

  comps <- data.frame(ComplexID=character(), ComplexName=character(),
                      subunitsUniprot_IDs=character(), subunitsNum=integer())
  data2 <- data[F,]
  for(h in 1:nrow(complexdb)) {
    subunitNames <- complexdb[h,"subunitsUniProt_IDs"]
    subunitname <- unique(unlist(strsplit(as.character(subunitNames),";")))
    #print(subunitname)
    for(s in subunitname) {
      pos <- grep(s, data$id, value=FALSE)
      # could have several isoforms, later to filter away the isoforms with fewer count number
      if (length(pos)>0) {
        data2 <- rbind(data2, data[pos, ])
        for (i in 1:length(pos)) {
          comps <- rbind(comps,complexdb[h,c(1,2,6,21)])
          #comps <- rbind(comps,complexdb[h,c(1,6,7,8)])
          # in the case of compleat database
        }
      }
    }
  }
  if (nrow(data2)==0) { stop("No protein complex found...") }
  data2 <- data2 %>% rowwise() %>% mutate(gene=getGeneName(description))
  comps3 <- cbind(comps,data2)
  # remove less quantified isoforms in one complex
  # remove complexes with too few subunits (number and percentage)
  comp_table <- comps3 %>% add_count(ComplexID,gene)
  cleariso <- FALSE
  if (sum(comp_table$n>1)) {
    print("Be cautious about the possible isoforms in the complex mapping...")
    print("Double check the potential_isoform_in_complex table in the working directory.")
    comp_table <- subset(comp_table, n>1)
    comp_table$n <- NULL
    # write.table(comp_table, paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"), dataname,
    #            "_potential_isoform_in_complex.txt"), sep="\t", row.names=FALSE, quote=FALSE)
    ms_filewrite(comp_table, paste0(dataname, "_potential_isoform_in_complex.txt"),
                 outdir=outdir, withdescription=F)
    cleariso <- TRUE
  }
  if (cutisoformbycount) {
    if (cleariso) { print("Would remove the isoforms with comparatively fewer count number.") }
    comps3 <- comps3 %>% left_join(data[,c("id","countNum")]) %>%
      group_by(ComplexID,gene) %>% top_n(1,countNum) %>%
      ungroup() %>% group_by(ComplexID) %>% mutate(subunitsIdentifiedNum=n()) %>%
      rowwise() %>% mutate(subunitsIdentifiedPerc=subunitsIdentifiedNum/subunitsNum) %>%
      filter(subunitsIdentifiedNum>=minsubunitsIdentified) %>% ungroup()
  }

  print(paste0("The number of identified multi-protein complexes with at least ",
               minsubunitsIdentified, " subunits is: ", length(unique(comps3$ComplexID))))

  return(comps3)
}
