function Fdeim = Fp_bad(F, V, vr, P)
    Ffull = F(V*vr);
    Fdeim = Ffull(P);
end