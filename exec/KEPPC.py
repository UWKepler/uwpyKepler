import sys

# koiTable_2012Feb26.txt 

infile = sys.argv[1]
for line in open(infile).readlines():
    KOI  = line[0:7].strip()
    KIC  = line[8:17].strip()
    Kp   = line[18:25].strip()
    t0   = line[26:36].strip()
    e_t  = line[37:45].strip()
    Per  = line[46:58].strip()
    e_P  = line[59:69].strip()
    Rad  = line[70:76].strip()
    a    = line[77:83].strip()
    Teq  = line[84:89].strip()
    Dur  = line[90:98].strip()
    Depth= line[99:105].strip()
    dR   = line[106:117].strip()
    e_dR = line[118:129].strip()
    rR   = line[130:138].strip()
    e_rR = line[139:147].strip()
    b    = line[148:156].strip()
    e_b  = line[157:165].strip()
    SNR  = line[166:172].strip()
    chi  = line[173:178].strip()
    Teff = line[179:184].strip()
    logg = line[185:190].strip()
    Rad  = line[191:197].strip()
    f_T  = line[198:200].strip()

    print ",".join((KOI, Dur, Depth, SNR, t0, e_t, Per, e_P, dR, e_dR, rR, e_rR, b, e_b, KIC, Teff, logg, Rad))

ingest = """
LOAD DATA INFILE '/data1/kepler/koiTable_2012Feb26.txt'
INTO TABLE KEPPC 
FIELDS TERMINATED BY ','
LINES TERMINATED BY '\n'
(KOI,Dur,Depth,SNR,Epoch,eEpoch,Period,ePeriod,aRstar,eaRstar,RpRstar,eRpRstar,b,eb,KID,Teff,logg,Rad);
"""
