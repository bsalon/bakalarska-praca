from arybo.lib import MBA

mba = MBA(1)

a = mba.var("a")
b = mba.var("b")
c = mba.var("c")
d = mba.var("d")

ab = a[0] * b[0]
ac = a[0] * c[0]
ad = a[0] * d[0]
bc = b[0] * c[0]
bd = b[0] * d[0]
cd = c[0] * d[0]
abc = a[0] * b[0] * c[0]
abd = a[0] * b[0] * d[0]
acd = a[0] * c[0] * d[0]
bcd = b[0] * c[0] * d[0]
abcd = a[0] * b[0] * c[0] * d[0]

polynom = mba.var("")

x = a + c + (a * b) + (b * d) + (c * d) + (a * b * d) + (b * c * d) + (a * b * c * d)
