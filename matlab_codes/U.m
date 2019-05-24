function val = U(tode)
    pp = spline(ptspan, U_, tode);
    val = pp(:,tode);
end