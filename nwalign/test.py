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

    s0 = "AGEBANAN"
    s1 = "ACEBAN"
    r = nw.global_align(s0, s1, gap_extend=-1, gap_open=-2, matrix="PAM250")
    assert r[0] == s0
    assert r[1] == s1 + "--"

    s0 = "ACEBAN"
    s1 = "AGEBANAN"
    r = nw.global_align(s0, s1, gap_extend=-1, gap_open=-2, matrix="PAM250")
    assert r[0] == s0 + "--"
    assert r[1] == s1

    s0 = "NABECA"
    s1 = "NANABEGA"
    r = nw.global_align(s0, s1, gap_extend=-1, gap_open=-2, matrix="PAM250")
    assert r[0] == "--" + s0
    assert r[1] == s1, r



    s0 = "AAAQDNVSLFYISAILNDMKEMPGIISRMPPLPVSINNDLASSLVTSATEPRN"
    s1 = "AAAEGDAWDHQPDAASNCLN"

    r = nw.global_align(s0, s1, gap_extend=-10, gap_open=-200, matrix='BLOSUM62')
    assert r[0] == s0
    assert r[1] == s1 + ("-" * (len(s0) - len(s1)))
    # using cached matrix.
    r2 = nw.global_align(s0, s1, gap_extend=-10, gap_open=-200, matrix='BLOSUM62')
    assert r == r2


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

if __name__ == "__main__":
    import nose
    nose.main()
