"""
Optimized kepler equation solver.
"""

@inline function kepler(e,M)
        M2 = M/2
        sinM2 = sin(M2)
        cosM2 = cos(M2)
        sinM22 = 2*sinM2
        one_e = 1 - e
        one_cosM = sinM22*sinM2
        sinM = sinM22*cosM2
        cosM = fma(-sinM22,sinM2,1.)
        a = e * sinM
        b = e * cosM
        one_b = fma(e,one_cosM,one_e) # 1 - e*cos(M) = (1-e) + e*(1-cos(M))

        # prediction step

        a2 = a^2

        c1 = fma(-1.565283593991395, b, 9.199102497914916)
        c2 = fma(c1, b, -18.479161685054063)
        c3 = fma(c2, b, 13.845342781130542)
        d1 = fma(-0.8547126454046926, b, 4.558693873335609)
        d2 = fma(d1,b,-5.068295197691118)
        d3 = fma(d2,a2,c3)
        r = a*d3

        c1 = fma(-0.26898397414324826, b, 2.1891776656780535)
        c2 = fma(c1, b, -5.571403408926362)
        c3 = fma(c2, b, 3.651209717391557)
        d1 = fma(-0.01700656596466221, b, 1.7445387962793053)
        d2 = fma(d1, b, -2.7065990285375796)
        q = fma(d2, a2, c3)

        c1 = fma(0.014326748319791131, b, 0.09959861662945457)
        c2 = fma(c1, b, 0.5164915422356884)
        c3 = fma(c2, b, -1.5282767152705534)
        d1 = fma(0.007184002827827581, b, 0.006245880306499962)
        d2 = fma(d1, b, 0.02923565692028258)
        d3 = fma(d2, a2, c3)
        p = a*d3


        # Solve cubic  (x-p)^3 + 3 * q (x-p) - 2 * r = 0
        aux = sqrt(fma(r,r,q*q*q)) + r  #   aux = sqrt(r^2 + q^3) + r
    #   s = sign(aux) * abs(aux)^0.33333333333333333
        s = sign(aux) * exp(0.3333333333333333333333*log(abs(aux))) # s = cbrt(aux)
        qs = q/s
        aux2 = fma(qs, qs, fma(s, s, q))
        aux3 = r / aux2
        x = fma(2.,aux3,p)


        # Correction step

        x2 = x*x

        # x - sin(x) Taylor polynomial, level 18

        v8 = fma(x2, 8.22063524662433e-18, -2.8114572543455206e-15)
        v7 = fma(x2, v8, 7.647163731819816e-13)
        v6 = fma(x2, v7, -1.6059043836821613e-10)
        v5 = fma(x2, v6, 2.505210838544172e-8)
        v4 = fma(x2, v5, -2.7557319223985893e-6)
        v3 = fma(x2, v4, 0.0001984126984126984)
        v2 = fma(x2, v3, -0.008333333333333333)
        v1 = fma(x2, v2, 0.1666666666666666666)
        bv1 = b*v1
        x_bsinx = x * fma(x2, bv1, one_b) # x-bsinx = (1-b)*x + b*(x - sinx)
        sinx = x * fma(-x2, v1, 1.) # sin(x) = x * (1 - x2 * v1)

        # 1 - cos(x) Taylor polynomial, level 18
        v8 = fma(x2, 1.5619206968586225e-16, -4.779477332387385e-14)
        v7 = fma(x2, v8, 1.1470745597729725e-11)
        v6 = fma(x2, v7, -2.08767569878681e-9)
        v5 = fma(x2, v6, 2.755731922398589e-7)
        v4 = fma(x2, v5, -2.48015873015873e-5)
        v3 = fma(x2, v4, 0.001388888888888889)
        v2 = fma(x2, v3, -0.041666666666666664)
        v1 = fma(x2, v2, 0.5)

        cosx = fma(-x2, v1, 1.) # cos(x) = 1. - x2 * v1
        one_cosx = v1 * x2
        one_bcosx = fma(b, one_cosx, one_b) # 1 - b*cos(x) = (1-b) + b*(1-cos(x))

        f = fma(-a, cosx, x_bsinx)
        df = fma(a, sinx, one_bcosx) # df = one_b + b * one_cosx + a * sinx
        df2 = fma(b, sinx, a * cosx)
        df3 = fma(b, cosx, - a * sinx)

        delta = - f / df
        eta = df2 / df
        nu = df3 / df

        aux2 = fma( delta*delta/3,-nu ,fma(eta,delta,2.) )
        aux3 = aux2 / (2*fma(eta,delta,1.))
        x = fma(delta, aux3, x)
        return M + x
end
