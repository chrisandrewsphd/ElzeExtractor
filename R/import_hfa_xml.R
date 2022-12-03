#' Import VF from HFA XML
#'
#' @param xmlfile Name of XML file to read
#' @param csvfilenamestub Base name of csv files. Protocol name will be appended.
#' @param idfunction Function to be applied to patient ID, for instance to correct errors or to de-identify (default: identity function)
#'
#' @return The protocol detected in the input file. (returned invisibly). Side effect of the function is to append the HVF values to the appropriate protocol file (if the file exists, or create a new file).
#' @export
#'
#' @examples
#' indir__ <- "S:/HFA Data/HFA2-Peds (Completed)"
#' outdir__ <- "C:/Users/UM_Ophthy_HSR/Downloads/hfa"
#' outname__ <- "pupil"
#'
#' xmlfiles__ <- dir(path=indir__, pattern="\\.xml$", full.names = TRUE, ignore.case = TRUE)
#' nfiles <- length(xmlfiles__)
#' cat("number of files to process:", nfiles, "\n")
#' protocols <- character(nfiles)

import.hfa.xml <- function(
		xmlfile,
		csvfilenamestub,
		idfunction = function(x) x)
{
	parsed = XML::xmlParse(xmlfile)
	patient = (XML::xmlChildren((XML::xmlChildren(parsed))$HFA_EXPORT))$PATIENT

	xg <- function(s, el=NULL, recursive=is.null(el))
	  XML::xmlElementsByTagName(if (is.null(el)) patient else el, s, recursive=recursive)
	xs <- function(...)
	{
		r = xg(...)
		ifelse(length(r)==0, "", XML::xmlValue((r[[1]])))
	}
	xn <- function(...) as.numeric(xs(...))

	protocol = sub(" Thr$", "", xs("DISPLAY_NAME"))		# measurement protocol without the trailing " Thr"
	if(grepl("^[0-9]+.*", protocol))	# sometimes "FT" is omitted
		protocol = paste0("FT-", protocol)

	if(!(protocol %in% c("SS-10-2", "SS-24-2", "SF-24-2", "SSW-24-2", "FP-24-2", "FT-24-2", "SS-30-2", "SF-30-2", "FP-30-2","FT-30-2", "FP-60-4")))
		#if(!(protocol %in% c("SS-10-2", "SSW-24-2")))
		warning(paste("read.hfa.xml: protocol", protocol, "not yet supported. Skipping file", xmlfile, "\n"))
	else {
		statpac = xn("STATPAC_STATUS") == 1	# if STATPAC not available: no global indices/TDs/PDs
		outfile = paste0(csvfilenamestub, "_", protocol, ifelse(statpac, "", "_no_statpac"), ".csv")
		ct <- function(..., comma=TRUE) cat(..., ifelse(comma, ",", "\n"), sep="", file=outfile, append=TRUE)
		repeated.measurements = !(substr(protocol, 1, 2) %in% c("SS", "SF"))		# SITA does not repeat locations
		nloc = xn("NUM_THRESHOLD_POINTS")
		c102tmp = c(1, 5, 7, 7, 9)
		##### generate hash table mapping VF locations to an index:
		loc.indices = as.data.frame(switch(paste(nloc),
																			 "54" = list(
																			 	x = c(-9, -3, 3, 9, -15, -9, -3, 3, 9, 15, -21, -15, -9, -3, 3, 9, 15, 21, -27, -21, -15, -9, -3, 3, 9, 15, 21, -27, -21, -15, -9, -3, 3, 9, 15, 21, -21, -15, -9, -3, 3, 9, 15, 21, -15, -9, -3, 3, 9, 15, -9, -3, 3, 9),
																			 	y = c(21, 21, 21, 21, 15, 15, 15, 15, 15, 15, 9, 9, 9, 9, 9, 9, 9, 9, 3, 3, 3, 3, 3, 3, 3, 3, 3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -9, -9, -9, -9, -9, -9, -9, -9, -15, -15, -15, -15, -15, -15, -21, -21, -21, -21)),	# 24-2
																			 "68" = list(x = unlist(sapply(-c(c102tmp, rev(c102tmp)), function(xleft) seq(xleft, -xleft, 2))),
																			 						y = unlist(mapply(rep.int, seq(9, -9, -2), c(c102tmp+1, rev(c102tmp+1))))), # 10-2
																			 "76" = list(
																			 	x = do.call("c", (sapply(c(-seq(9, 27, 6), -27, -27, seq(-27,-9,6)), function(xleft) seq(xleft, -xleft, 6)))),
																			 	y = c(rep(27,4), rep(21,6), rep(15,8), rep(9,10), rep(3,10), rep(-3,10), rep(-9,10), rep(-15,8), rep(-21,6), rep(-27,4))),		# 30-2
																			 "60" = list(
																			 	x = c(seq(-30,30,12), seq(-54,54,12), rep(c(seq(-54,-30,12), seq(30, 54, 12)), 4), rep(seq(-42,42,12), 2), seq(-18, 18, 12)),
																			 	y =  c(rep(42,6), rep(30,10), rep(18,6), rep(6,6), rep(-6,6), rep(-18,6), rep(-30,8), rep(-42,8), rep(-54, 4)))))		# 60-4
		vfloc2ind.ht = new.env()
		for(i in 1:nrow(loc.indices))
			vfloc2ind.ht[[paste(loc.indices[i,], collapse=",")]] <- i

		exclude.blindspot <- function(v)
			if(nloc == 54 || nloc == 76)
			{
				if(nloc == 54)
					v[-c(26,35)]
				else
					v[-c(36, 46)]
			} else v

		numstring <- function(prefix="s", blindspot=substr(prefix, 1, 1) != 's')
		{
			v = paste0(prefix, 1:nloc)
			r = if(blindspot) exclude.blindspot(v) else v
			paste(r, collapse=",")
		}

		if(!file.exists(outfile))
		{
			# generate file and write header:
			header = paste0(
				#"id,age,righteye,timeoftest,duration,centralval,centralprob,pupil,sphere,cylinder,axis,distsphere,distcylinder,distaxis,acuity,",
				"id,age,righteye,dateoftest,duration,centralval,centralprob,pupil,sphere,cylinder,axis,distsphere,distcylinder,distaxis,acuity,",
				ifelse(repeated.measurements, "falsenegnum,falsenegdenom,falseposnum,falseposdenom", "falsenegrate,falseposrate"),
				",malfixnum,malfixdenom,blindspotx,blindspoty,ght,vfi,md,mdprob,psd,psdprob,",
				numstring(), ",",
				ifelse(repeated.measurements, paste0(numstring("sa"), ","), ""),
				ifelse(
					statpac,
					paste(numstring("td"), numstring("pd"), numstring("tdp"), numstring("pdp"), sep=","),
					numstring("dd")),	# dd: defect depth (if no STATPAC available, such as for stimulus size V)
				"\n")
			cat(header, file=outfile)
		}

		### write values to outfile:
		ptid = idfunction(xs("PATIENT_ID", recursive=FALSE))
		ct(ptid)

		dob <- as.Date(xs("BIRTH_DATE", recursive=FALSE), format = "%Y-%m-%d")
		study = XML::xmlElementsByTagName(patient, "STUDY")[[1]]
		vid <- as.Date(xs("VISIT_DATE", study), format = "%Y-%m-%d")
		age = round(as.numeric(vid-dob)/365.25, 2)
		ct(age)

		series = XML::xmlElementsByTagName(study, "SERIES")[[1]]
		righteye = xs("SITE", series)
		oscoef = ifelse(righteye=="0", -1, 1)
		ct(righteye)

		# td = xs("SERIES_DATE_TIME", series)
		# timeoftest = paste(
		# 	as.numeric(substr(td, 3, 4)),
		# 	substr(td, 6, 7),
		# 	ifelse(nchar(td) < 19, paste("0", substr(td, 9,9), sep=""), substr(td, 9,10)),
		# 	ifelse(nchar(td) < 19, substr(td, 11,12), substr(td, 12,13)),
		# 	ifelse(nchar(td) < 19, substr(td, 14,15), substr(td, 15,16)),
		# 	sep="")
		# ct(timeoftest)
		ct(as.character(vid))

		exam = xg("FIELD_EXAM", series)[[1]]
		vftest = xg("STATIC_TEST", exam)[[1]]
		dur = xs("EXAM_DURATION", vftest)
		duration = 60*as.numeric(substr(dur, 4, 5)) + as.numeric(substr(dur, 7, 8))
		ct(duration)

		ct(xn("FOVEAL_THRESHOLD", vftest))
		ct(xn("CENTRAL_REF_LEVEL", vftest))

		ct(xn("PUPIL_DIAMETER", exam)) # new search position

		triallensxml = xg("TRIAL_RX", exam)[[1]]
		ct(xn("SPHERE", triallensxml))
		ct(xn("CYLINDER", triallensxml))
		ct(xn("AXIS", triallensxml))
		distlensxml = xg("DISTANCE_RX", exam)[[1]]
		ct(xn("SPHERE", distlensxml))
		ct(xn("CYLINDER", distlensxml))
		ct(xn("AXIS", distlensxml))

		ct(xn("VA_STRING", exam))

		print.rel.index <- function(s)
		{
			ind = xg(s, vftest)[[1]]
			num = xn("ERRORS", ind)
			denom = xn("TRIALS", ind)
			ct(num, ",", denom)
		}

		if(repeated.measurements)
			lapply(c("FALSE_NEGATIVES", "FALSE_POSITIVES"), print.rel.index)
		else {
			print.rate <- function(s)
			{
				r = xn(s, vftest)/100
				ct(ifelse(r==-1, NA, r))
			}
			lapply(c("FALSE_NEGATIVE_PERCENT", "FALSE_POSITIVE_PERCENT"), print.rate)
		}
		print.rel.index("FIXATION_CHECK")
		ct(xn("BLIND_SPOT_X", vftest))
		ct(xn("BLIND_SPOT_Y", vftest))

		thresh = xg("THRESHOLD_TEST", vftest)[[1]]

		# global indices:
		if(statpac)
		{
			stpc = xg("STATPAC", thresh)[[1]]
			ct(xn("GHT", stpc))
			glind = xg("GLOBAL_INDICES", stpc)[[1]]
			print.ind <- function(s)
				ct(xn(s, glind))
			lapply(c("VFI", "MD", "MD_PROBABILITY", "PSD", "PSD_PROBABILITY"), print.ind)
		} else
			ct(paste(rep(NA, 6), collapse=","))

		# sensitivities:
		sensxml = xg("THRESHOLD_XY_LOCATION", xg("THRESHOLD_SITE_LIST", thresh)[[1]])
		sens = rep(NA, nloc)
		sensrep = sens
		assign.sensitivities <- function(el)
		{
			g <- function(s)
				as.numeric(XML::xmlValue(XML::xmlElementsByTagName(el, s)[[1]]))
			get.value <- function(rep=1)
			{
				r = g(paste0("RESULT_", rep))	# result=2 means "<0", decoded as -2 here
				ifelse(r==2, -2, g(paste0("THRESHOLD_", rep)))
			}
			loc = paste(oscoef*g("X"), g("Y"), sep=",")
			sens[ vfloc2ind.ht[[loc]] ] <<- get.value(1)
			if(repeated.measurements && length(XML::xmlChildren(el))>4)
				sensrep[ vfloc2ind.ht[[loc]] ] <<- get.value(2)
		}
		lapply(sensxml, assign.sensitivities)
		ct(paste(sens, collapse=","))
		if(repeated.measurements)
			ct(paste(sensrep, collapse=","))

		# TDs/PDs or defect depths:
		print.normative.values <- function(s="TOTAL_DEVIATION", ...)
		{
			v = rep(NA, nloc)
			listname = switch(s,
												"DEFECT_DEPTH" = "DEFECT_DEPTH_SITE_LIST",
												"TOTAL_DEVIATION_PROBABILITY" = "TOTAL_DEVIATION_PROBABILITY_LIST",
												"PATTERN_DEVIATION_PROBABILITY" = "PATTERN_DEVIATION_PROBABILITY_LIST",
												paste0(s, "_VALUE_LIST"))

			locname = ifelse(grepl(".*PROBABILITY$", s), s, paste0(s, "_VALUE"))

			vallist = XML::xmlChildren(xg(listname, thresh, recursive=TRUE)[[1]])
			assign.val <- function(el)
			{
				g <- function(s)
					as.numeric(XML::xmlValue(XML::xmlElementsByTagName(el, s)[[1]]))
				loc = paste(oscoef*g("X"), g("Y"), sep=",")
				v[ vfloc2ind.ht[[loc]] ] <<- g(locname)
			}
			lapply(vallist, assign.val)
			r = exclude.blindspot(v)
			ct(paste(r, collapse=","), ...)
		}

		if(statpac)
		{
			print.normative.values("TOTAL_DEVIATION")
			print.normative.values("PATTERN_DEVIATION")
			print.normative.values("TOTAL_DEVIATION_PROBABILITY")
			print.normative.values("PATTERN_DEVIATION_PROBABILITY", comma=FALSE)
		} else
			print.normative.values("DEFECT_DEPTH", comma=FALSE)
	}
	#return(invisible(NULL))
	return(invisible(protocol))
}


#########
# usage example (with your attached XML files);
# place this script into the same directory as the XML files
# or set the directory # e.g.
# setwd("../Input")
# can specify output also
#########

# # default
# indir__ <- getwd()
# outdir__ <- indir__
#
# # example: get all xml files in current directory:
# xmlfiles__ <- dir(path=indir__, pattern="\\.xml$")
# lapply(
# 	xmlfiles__,
# 	function(filename) import.hfa.xml(xmlfile = filename, csvfilenamestub = "hfa_measurements"))
#

#
#
#
#
#
# outbase__ <- sprintf("%s/%s", outdir__, outname__)
# system.time(
# 	for (i in 1:100) {#seq.int(nfiles)) {
# 		if (!(i %% 1000)) cat("file ", i, "\n")
# 		protocols[i] <- import.hfa.xml(xmlfile = xmlfiles__[i], csvfilenamestub = outbase__)
# 	}
# )
# head(protocols)
# t(t(table(protocols)))
#
#
#
