# Estimators.

function estimator1(R::Float64, beta::Float64, int::Integrator,
                    s::State{P,D,A,N}) where {P,D,A,N}
    function estim()
        result = 2.0 / R # 1 / nm

        dxdr = 0.5 .* normalize(com_disp(CN1, CN2, CJ, s))

        for a in 1:A
            for d in 1:D
                for j in 1:P
                    result += dxdr[d] * get_force(int)[j, d, a, CN1] * beta
                    result -= dxdr[d] * get_force(int)[j, d, a, CN2] * beta
                end
            end
        end

        result
    end
end

function estimator2(R::Float64, beta::Float64, int::Integrator,
                    s::State{P,D,A,N}) where {P,D,A,N}
    function estim()
        result = 2.0 / R # 1 / nm

        dxdr = 0.5 .* normalize(com_disp(CN1, CN2, CJ, s))

        for a in 1:A
            for d in 1:D
                result += dxdr[d] * get_force(int)[CJ, d, a, CN1] * beta
                result -= dxdr[d] * get_force(int)[CJ, d, a, CN2] * beta

                diff1 = (2 * s.qs[CJ, d, a, CN1]
                         - s.qs[bead_next(P, CJ), d, a, CN1]
                         - s.qs[bead_prev(P, CJ), d, a, CN1]) # nm
                result -= dxdr[d] * diff1 * MS[a] * P / (hbar^2 * beta)
                diff2 = (2 * s.qs[CJ, d, a, CN2]
                         - s.qs[bead_next(P, CJ), d, a, CN2]
                         - s.qs[bead_prev(P, CJ), d, a, CN2]) # nm
                result += dxdr[d] * diff2 * MS[a] * P / (hbar^2 * beta)
            end
        end

        result
    end
end
