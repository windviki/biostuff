import nwalign as nw
print nw

def test_nw():

    r = nw.global_align("CEOLECANTH", "PELICAN")
    assert r ==  ('CEOLECANTH', 'PE-LICAN--'), r
    r = nw.global_align("PELICAN", "CEOLECANTH")
    assert r ==  ('PE-LICAN--', 'CEOLECANTH'), r

def test_pam():
    r = nw.global_align("CEELECANTHHH", "PELICAN", gap_open=-2, gap_extend=-1, matrix='PAM250')
    assert r ==  ('CEELECANTHHH', '-PELICAN----'), r
    r = nw.global_align("PELICAN", "CEELECANTHHH", gap_open=-2, gap_extend=-1, matrix='PAM250')
    assert r ==  ('-PELICAN----', 'CEELECANTHHH'), r

def test_score():
    s0 = "AGEBANAN"
    s1 = "ACEBAN"
    r = nw.global_align(s0, s1, gap_extend=-1, gap_open=-2, matrix="PAM250")
    s = nw.score_alignment(r[0], r[1], gap_open=-2, gap_extend=-1, matrix="PAM250")
    assert s == 7, s

def test_gap_open():

    s0 = "ACEBANAN"
    s1 = "ACEBAN"
    r = nw.global_align(s0, s1, gap_extend=-1, gap_open=-2, matrix="PAM250")
    print r[0]
    print r[1]
    assert r[0] == s0, r[0]
    assert r[1] == s1 + "--", r

    s0 = "ACEBAN"
    s1 = "ASEBANAN"
    r = nw.global_align(s0, s1, gap_extend=-1, gap_open=-2, matrix="PAM250")
    assert r[0] == s0 + "--", r
    assert r[1] == s1, r


    s0 = "WWWWQDNVSLFYISAILNDMKEMPGIISRMPPLPVSINNDLASSLVTSATEPRN"
    s1 = "WWWWEGDAWDHQPDAASNCLN"

    s0 = "WWWWQD"
    s1 = "WWWWEDDA"

    r = nw.global_align(s0, s1, gap_extend=-10, gap_open=-20, matrix='BLOSUM62')
    print r[0]
    print r[1]
    assert r[0] == s0 + "--"

    assert r[1] == s1 
    # using cached matrix.
    r2 = nw.global_align(s0, s1, gap_extend=-10, gap_open=-20, matrix='BLOSUM62')
    assert r == r2


def test_marcin():
    r = nw.global_align('CPEL', 'PREK', gap_open=-6,
gap_extend=-2, matrix='BLOSUM62')
    assert r[0] == "CPEL"
    assert r[1] == "PREK"

def test_raises():
    a="TTAAT"
    b="TT"

    al = "TTAAT"
    bl = "TT---"

    r = nw.global_align(a, b, gap_open=-2, gap_extend=-1)
    assert r[0] == al, r[0]
    assert r[1] == bl, r[1]

    from nose.tools import assert_raises
    assert_raises(AssertionError, nw.global_align, a, b, gap_open=2)

def test_bad_char():
    a="tTAAT&"
    b="TT"
    r2 = nw.global_align(a, b, gap_extend=-10, gap_open=-20, matrix='BLOSUM62')

def test_no_string():
    a = ""
    b = "TT"
    r2 = nw.global_align("", "TT", gap_extend=-10, gap_open=-20, matrix='BLOSUM62')
    r2 = nw.global_align("TT", "", gap_extend=-10, gap_open=-20, matrix='BLOSUM62')
    r2 = nw.global_align("", "", gap_extend=-10, gap_open=-20, matrix='BLOSUM62')
    r2 = nw.global_align("", "TT", gap_extend=-10, gap_open=-20)
    r2 = nw.global_align("TT", "", gap_extend=-10, gap_open=-20)
    r2 = nw.global_align("", "", gap_extend=-10, gap_open=-20)

if __name__ == "__main__":
    import nose
    nose.main()
