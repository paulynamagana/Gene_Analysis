library(SRAdb)
sqlfile <- getSRAdbFile()
sra_con <- dbConnect(SQLite(),sqlfile)
rs = getSRAfile("SRR000657", sra_con, fileType = 'sra' )
