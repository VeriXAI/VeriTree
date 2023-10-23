@commutative function minkowski_sum(B1::OpenHyperrectangle,
                                    B2::AbstractHyperrectangle)
    return OpenHyperrectangle(B1.lo + low(B2), B1.hi + high(B2),
                              B1.lo_open, B1.hi_open)
end
