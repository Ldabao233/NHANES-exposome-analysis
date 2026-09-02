#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(haven)
  library(tibble)
})

xpt_dir <- Sys.getenv("XPT_DIR", ".")
out_dir <- Sys.getenv("RDS_DIR", ".")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_merge <- function(files) {
  x <- rbindlist(lapply(file.path(xpt_dir, files), read_xpt), fill = TRUE)
  as_tibble(x[, lapply(.SD, function(v) na.omit(v)[1]), by = SEQN])
}

build_module <- function(files, map, output, urine = FALSE, scale_factor = 100,
                         multiply_1000 = character()) {
  raw <- read_merge(files)
  if (urine) raw <- raw |> select(-any_of("URXUCR")) |> left_join(creatinine, by = "SEQN")
  result <- tibble(SEQN = raw$SEQN)

  for (target in names(map)) {
    source <- intersect(strsplit(map[[target]], "|", fixed = TRUE)[[1]], names(raw))
    stopifnot(length(source) > 0)
    result[[target]] <- Reduce(coalesce, raw[source])
  }

  result[multiply_1000] <- lapply(result[multiply_1000], function(x) x * 1000)

  if (urine) {
    analytes <- names(map)[!endsWith(names(map), "LC")]
    result[analytes] <- lapply(result[analytes], function(x) x / raw$URXUCR * scale_factor)
    result$creatinine <- raw$URXUCR
  }

  saveRDS(result, file.path(out_dir, output))
}

creatinine <- read_merge(c(
  "ALB_CR_H_2013-2014.xpt",
  "ALB_CR_I_2015-2016.xpt"
)) |>
  select(SEQN, URXUCR)

# Serum flame retardants
years <- c("2013-2014", "2015-2016")
pool_ids <- bind_rows(
  read_xpt(file.path(xpt_dir, "POOLTF_H.xpt")) |> mutate(year = years[1]),
  read_xpt(file.path(xpt_dir, "POOLTF_I.xpt")) |> mutate(year = years[2])
) |>
  select(SEQN, SAMPLEID, year)
bfr_pool <- bind_rows(
  read_xpt(file.path(xpt_dir, "BFRPOL_H.xpt")) |> mutate(year = years[1]),
  read_xpt(file.path(xpt_dir, "BFRPOL_I.xpt")) |> mutate(year = years[2])
) |>
  left_join(pool_ids, by = c("SAMPLEID", "year"))

BFR_serum <- bfr_pool |>
  transmute(
    SEQN,
    PBB153 = LBCBB1LA, PBB153LC = LBDBB1LC,
    PBDE28 = LBCBR2LA, PBDE28LC = LBDBR2LC,
    PBDE47 = LBCBR3LA, PBDE47LC = LBDBR3LC,
    PBDE85 = LBCBR4LA, PBDE85LC = LBDBR4LC,
    PBDE99 = LBCBR5LA, PBDE99LC = LBDBR5LC,
    PBDE100 = LBCBR6LA, PBDE100LC = LBDBR6LC,
    PBDE153 = LBCBR7LA, PBDE153LC = LBDBR7LC,
    PBDE154 = LBCBR8LA, PBDE154LC = LBDBR8LC,
    PBDE183 = LBCBR9LA, PBDE183LC = LBDBR9LC,
    PBDE209 = LBCBR11L, PBDE209LC = LBDBR11C
  )
saveRDS(BFR_serum, file.path(out_dir, "BFR_serum.rds"))

build_module(
  c("PBCD_H.xpt", "IHGEM_H.xpt", "PBCD_I.xpt", "IHGEM_I.xpt"),
  c(
    Cd = "LBDBCDSI", CdLC = "LBDBCDLC",
    Pb = "LBDBPBSI", PbLC = "LBDBPBLC",
    THg = "LBDTHGSI", THgLC = "LBDTHGLC",
    IHg = "LBDIHGSI", IHgLC = "LBDIHGLC",
    Se = "LBDBSESI", SeLC = "LBDBSELC",
    Mn = "LBDBMNSI", MnLC = "LBDBMNLC"
  ),
  "HM_serum.rds"
)

build_module(
  c("PFC_Serum_2013-2014.xpt", "PFC_Serum_2015-2016.xpt"),
  c(
    n_PFOA = "LBXNFOA", n_PFOALC = "LBDNFOAL",
    n_PFOS = "LBXNFOS", n_PFOSLC = "LBDNFOSL",
    Sm_PFOS = "LBXMFOS", Sm_PFOSLC = "LBDMFOSL",
    PFDA = "LBXPFDE", PFDALC = "LBDPFDEL",
    PFHS = "LBXPFHS", PFHSLC = "LBDPFHSL",
    PFHpA = "LBXPFHP", PFHpALC = "LBDPFHPL",
    PFNA = "LBXPFNA", PFNALC = "LBDPFNAL",
    PFUnDA = "LBXPFUA", PFUnDALC = "LBDPFUAL",
    MeFOSAA = "LBXMPAH", MeFOSAALC = "LBDMPAHL"
  ),
  "PFC_serum.rds"
)

build_module(
  c("VOCWB_H.xpt", "VOCWB_I.xpt"),
  c(
    X25_DMF = "LBX2DF", X25_DMFLC = "LBD2DFLC",
    Benzene = "LBXVBZ", BenzeneLC = "LBDVBZLC",
    X14_DCB = "LBXVDB", X14_DCBLC = "LBDVDBLC",
    Ethylbenzene = "LBXVEB", EthylbenzeneLC = "LBDVEBLC",
    Furan = "LBXVFN", FuranLC = "LBDVFNLC",
    o_Xylene = "LBXVOX", o_XyleneLC = "LBDVOXLC",
    m_p_Xylene = "LBXVXY", m_p_XyleneLC = "LBDVXYLC",
    BDCM = "LBXVBM", BDCMLC = "LBDVBMLC",
    Chloroform = "LBXVCF", ChloroformLC = "LBDVCFLC",
    DBCM = "LBXVCM", DBCMLC = "LBDVCMLC",
    Toluene = "LBXVTO", TolueneLC = "LBDVTOLC"
  ),
  "VOC_serum.rds"
)

build_module(
  c("ETHOX_H.xpt", "FORMAL_H.xpt", "ETHOX_I.xpt", "FORMAL_I.xpt"),
  c(EO = "LBXEOA", EOLC = "LBDEOALC", FA = "LBXFOR", FALC = "LBDFORLC"),
  "gas_serum.rds"
)

build_module(
  c("SSFLRT_H.xpt", "SSFR_I.xpt"),
  c(
    DPHP = "SSDPHP|URXDPHP", DPHPLC = "SSDPHPL|URDDPHLC",
    BDCPP = "SSBDCPP|URXBDCP", BDCPPLC = "SSBDCPPL|URDBDCLC",
    BCPP = "SSBCPP|URXBCPP", BCPPLC = "SSBCPPL|URDBCPLC",
    BCETP = "SSBCEP|URXBCEP", BCETPLC = "SSBCEPL|URDCEPLC",
    DBUP = "SSDBUP|URXDBUP", DBUPLC = "SSDBUPL|URDDUPLC",
    DPCP = "SSDPCP", DPCPLC = "SSDPCPL",
    IPPPP = "SSIPPP", IPPPPLC = "SSIPPPL",
    TBPPP = "SSBPPP", TBPPPLC = "SSBPPPL"
  ),
  "BFR_urine.rds", TRUE
)

build_module(
  c("Phthalates_2013-2014.xpt", "Phthalates_2013-2014-2.xpt",
    "Phthalates_2015-2016.xpt"),
  c(
    MHBP = "SSURMHBP|URXMHBP", MHBPLC = "SDUMHBPL|URDMHBLC",
    HIBP = "SSURHIBP|URXHIBP", HIBPLC = "SDUHIBPL|URDHIBLC",
    MnBP = "URXMBP", MnBPLC = "URDMBPLC",
    MEP = "URXMEP", MEPLC = "URDMEPLC",
    MEHP = "URXMHP", MEHPLC = "URDMHPLC",
    MINP = "URXMNP", MINPLC = "URDMNPLC",
    MBzP = "URXMZP", MBzPLC = "URDMZPLC",
    MCPP = "URXMC1", MCPPLC = "URDMC1LC",
    MEHHP = "URXMHH", MEHHPLC = "URDMHHLC",
    MEOHP = "URXMOH", MEOHPLC = "URDMOHLC",
    MiBP = "URXMIB", MiBPLC = "URDMIBLC",
    MECPP = "URXECP", MECPPLC = "URDECPLC",
    MCOCH = "URXMCOH", MCOCHLC = "URDMCOLC",
    MCNP = "URXCNP", MCNPLC = "URDCNPLC",
    MCOP = "URXCOP", MCOPLC = "URDCOPLC",
    MHNCH = "URXMHNC", MHNCHLC = "URDMCHLC"
  ),
  "Phthalates_urine.rds", TRUE
)

build_module(
  c("DEET_2013-2014.xpt", "UPHOPM_2013-2014.xpt",
    "DEET_2015-2016.xpt", "UPHOPM_2015-2016.xpt",
    "OPI_2015-2016.xpt", "NEON_2015-2016.xpt"),
  c(
    DEA = "URXDEA", DEALC = "URDDEALC",
    DMP = "URXOP1", DMPLC = "URDOP1LC",
    DEP = "URXOP2", DEPLC = "URDOP2LC",
    DMTP = "URXOP3", DMTPLC = "URDOP3LC",
    DETP = "URXOP4", DETPLC = "URDOP4LC",
    DMDTP = "URXOP5", DMDTPLC = "URDOP5LC",
    X24D = "URX24D", X24DLC = "URD24DLC",
    X4FPBA = "URX4FP", X4FPBALC = "URD4FPLC",
    TCPY = "URXCPM", TCPYLC = "URDCPMLC",
    MDA = "URXMAL", MDALC = "URDMALLC",
    X3PBA = "URXOPM", X3PBALC = "URDOPMLC",
    OXY = "URXOXY", OXYLC = "URDOXYLC",
    PNP = "URXPAR", PNPLC = "URDPARLC",
    DCVDMPCA = "URXTCC", DCVDMPCALC = "URDTCCLC",
    X5OHIMI = "SSOHIM", X5OHIMILC = "SSOHIMLC",
    NDMA = "SSAND", NDMALC = "SSANDLC"
  ),
  "pesticides_urine.rds", TRUE
)

build_module(
  c("EP_PP_2013-2014.xpt", "EP_PP_2015-2016.xpt"),
  c(
    BPA = "URXBPH", BPALC = "URDBPHLC",
    BP3 = "URXBP3", BP3LC = "URDBP3LC",
    BuP = "URXBUP", BuPLC = "URDBUPLC",
    EtP = "URXEPB", EtPLC = "URDEPBLC",
    MeP = "URXMPB", MePLC = "URDMPBLC",
    PrP = "URXPPB", PrPLC = "URDPPBLC",
    BPF = "URXBPF", BPFLC = "URDBPFLC",
    BPS = "URXBPS", BPSLC = "URDBPSLC",
    TCC = "URXTLC", TCCLC = "URDTLCLC",
    X25_DCP = "URX14D", X25_DCPLC = "URD14DLC",
    X24_DCP = "URXDCB", X24_DCPLC = "URDDCBLC"
  ),
  "PCCP_urine.rds", TRUE
)

build_module(
  c("SSGLYP_2013-2014.xpt", "SSGLYP_2015-2016.xpt"),
  c(GLYP = "SSGLYP", GLYPLC = "SSGLYPL"),
  "GYLP_urine.rds", TRUE
)

build_module(
  c("PFC_urine_2013-2014.xpt"),
  c(PFBA = "SSPFBA", PFBALC = "SSPFBAL",
    PFHxA = "SSPFHA", PFHxALC = "SSPFHAL"),
  "PFC_urine.rds", TRUE
)

build_module(
  c("AA_2013-2014.xpt", "HCAA_2013-2014.xpt", "UADM_2015-2016.xpt"),
  c(
    X1_AN = "URX1NP", X1_ANLC = "URD1NPLC",
    X2_AN = "URX2NP", X2_ANLC = "URD2NPLC",
    X4_ABP = "URX4BP", X4_ABPLC = "URD4BPLC",
    o_ANS = "URXANS", o_ANSLC = "URDANSLC",
    X26_DMA = "URXDMN", X26_DMALC = "URDDMNLC",
    o_TLD = "URXOTD", o_TLDLC = "URDOTDLC",
    X4_MDA = "URX4MDA", X4_MDALC = "URD4MALC",
    PPDA = "URXPPDA", PPDALC = "URDPDALC",
    A_a_C = "URXAAC", A_a_CLC = "URDAACLC",
    Harman = "URXHM", HarmanLC = "URDHMLC",
    MeA_a_C = "URXMAAC", MeA_a_CLC = "URDMACLC",
    Norharman = "URXNHM", NorharmanLC = "URDNHMLC",
    PhIP = "URXPHIP", PhIPLC = "URDPHPLC"
  ),
  "AA_urine.rds", TRUE, 10, c("X4_MDA", "PPDA")
)

build_module(
  c("PAH_2013-2014.xpt", "PAH_2015-2016.xpt"),
  c(
    X1_NP = "URXP01", X1_NPLC = "URDP01LC",
    X2_NP = "URXP02", X2_NPLC = "URDP02LC",
    X3_FL = "URXP03", X3_FLLC = "URDP03LC",
    X2_FL = "URXP04", X2_FLLC = "URDP04LC",
    X1_PH = "URXP06", X1_PHLC = "URDP06LC",
    X1_PY = "URXP10", X1_PYLC = "URDP10LC",
    X23_OH_PH = "URXP25", X23_OH_PHLC = "URDP25LC"
  ),
  "PAH_urine.rds", TRUE, 10
)

build_module(
  c("UHM_H.xpt", "UAS_H.xpt", "UHG_H.xpt", "UTAS_H.xpt",
    "UHM_I.xpt", "UAS_I.xpt", "UHG_I.xpt", "UTAS_I.xpt"),
  c(
    Hg = "URXUHG", HgLC = "URDUHGLC",
    Ba = "URXUBA", BaLC = "URDUBALC",
    Cd = "URXUCD", CdLC = "URDUCDLC",
    Co = "URXUCO", CoLC = "URDUCOLC",
    Cs = "URXUCS", CsLC = "URDUCSLC",
    Mo = "URXUMO", MoLC = "URDUMOLC",
    Pb = "URXUPB", PbLC = "URDUPBLC",
    Sb = "URXUSB", SbLC = "URDUSBLC",
    Tl = "URXUTL", TlLC = "URDUTLLC",
    W = "URXUTU", WLC = "URDUTULC",
    U = "URXUUR", ULC = "URDUURLC",
    Mn = "URXUMN", MnLC = "URDUMNLC",
    Sn = "URXUSN", SnLC = "URDUSNLC",
    Sr = "URXUSR", SrLC = "URDUSRLC",
    Total_As = "URXUAS", Total_AsLC = "URDUASLC",
    As_III = "URXUAS3", As_IIILC = "URDUA3LC",
    AsC = "URXUAC", AsCLC = "URDUACLC",
    AsB = "URXUAB", AsBLC = "URDUABLC",
    DMA = "URXUDMA", DMALC = "URDUDALC",
    MMA = "URXUMMA", MMALC = "URDUMMAL"
  ),
  "HM_urine.rds", TRUE
)

build_module(
  c("UVOC_H.xpt", "UVOC_I.xpt"),
  c(
    X2_MHA = "URX2MH", X2_MHALC = "URD2MHLC",
    X34_MHA = "URX34M", X34_MHALC = "URD34MLC",
    CEMA_AAM = "URXAAM", CEMA_AAMLC = "URDAAMLC",
    CMAC = "URXAMC", CMACLC = "URDAMCLC",
    ATCA = "URXATC", ATCALC = "URDATCLC",
    BMA = "URXBMA", BMALC = "URDBMALC",
    NPMA = "URXBPM", NPMALC = "URDBPMLC",
    CEMA_CEM = "URXCEM", CEMA_CEMLC = "URDCEMLC",
    CYHA_Cys = "URXCYHA", CYHA_CysLC = "URDCYALC",
    CEMA_CYM = "URXCYM", CEMA_CYMLC = "URDCYMLC",
    DHB_Cys = "URXDHB", DHB_CysLC = "URDDHBLC",
    GAM_Cys = "URXGAM", GAM_CysLC = "URDGAMLC",
    HEMA = "URXHEM", HEMALC = "URDHEMLC",
    X2_HPMA = "URXHP2", X2_HPMALC = "URDHP2LC",
    X3_HPMA = "URXHPM", X3_HPMALC = "URDHPMLC",
    IPM1_Cys = "URXIPM1", IPM1_CysLC = "URDPM1LC",
    IPM3_Cys = "URXIPM3", IPM3_CysLC = "URDPM3LC",
    MA = "URXMAD", MALC = "URDMADLC",
    MB3_Cys = "URXMB3", MB3_CysLC = "URDMB3LC",
    PHE_Cys = "URXPHE", PHE_CysLC = "URDPHELC",
    PGA = "URXPHG", PGALC = "URDPHGLC",
    PMA = "URXPMA", PMALC = "URDPMALC",
    PMM_Cys = "URXPMM", PMM_CysLC = "URDPMMLC",
    TTCA = "URXTTC", TTCALC = "URDTTCLC"
  ),
  "VOC_urine.rds", TRUE
)

build_module(
  c("PERNT_H.xpt", "PERNT_I.xpt"),
  c(
    Perchlorate = "URXUP8", PerchlorateLC = "URDUP8LC",
    Nitrate = "URXNO3", NitrateLC = "URDNO3LC",
    Thiocyanate = "URXSCN", ThiocyanateLC = "URDSCNLC"
  ),
  "PERNT_urine.rds", TRUE
)
