rm(list = ls())

cmodel <- "FATES"

origin.dir <- "/home/femeunier/Documents/projects/Yoko.regrowth/outputs/raw.FATES/"
dest.dir <- "/home/femeunier/Documents/projects/Yoko.regrowth/outputs/MIP/FATES/"


files <- list.files(origin.dir,
                    pattern = ".*cVegAbov.*.nc$",
                    recursive = TRUE,full.names = TRUE)

for (ifile in seq(1,length(files))){

  cfile <- tools::file_path_sans_ext(basename(files[ifile]))
  cfile_split1 <- strsplit(tools::file_path_sans_ext((cfile)),"\\_")[[1]]
  cfile_split1[length(cfile_split1)] <- "cAGB"


  dest.name <- paste0(paste0(cfile_split1, collapse = "_"),
                      ".nc")

  print(dest.name)

  system2("nccopy",
          c("-k nc4",files[ifile],
            file.path(dest.dir,dest.name)))

}
