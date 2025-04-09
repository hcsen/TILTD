function [perturbed, optim_archive, m_archive, violation_archive] = basinhop(k, optimised, sigmas, MBH_tail, ...
        MBH_theta, rho_hop, t0Hop, dtHop, dtIndex, lb, ub, ind_vinfi, ind_vinff, s_v, N_flybys, ...
        phaseSizes, objInd, whichThrust, consts, options, A, b, Aeq, beq, Np, optim_archive, m_archive, violation_archive)

        probSize = length(optimised);

        % Perturb all the decision variables
        pm = fix(rand(1, probSize) + 0.5);
        pm(~pm) = -1;
        perturbed = optimised + gprnd(MBH_tail, sigmas, MBH_theta * ones(1, probSize)) .* pm;

        % Hop to a new random guess if needed
        if rand < rho_hop
            disp(k);
            perturbed(1) = perturbed(1) + 2 * t0Hop * rand - t0Hop; % Hop the launch epoch
            for j = 1:Np
                perturbed(dtIndex(j)) = perturbed(dtIndex(j)) + 2 * dtHop * rand - dtHop; % Hop each phase tof
            end
        end

        % Ensure decision variables are within bounds
        for j = 1:probSize
            if perturbed(j) < lb(j)
                perturbed(j) = lb(j);
                % Adjust vinfi if it equals vinff
                if N_flybys > 0
                    for vind = 1:N_flybys
                        if norm(perturbed(ind_vinfi((1:3) + (vind - 1) * 3))) == norm(perturbed(ind_vinff((1:3) + (vind - 1) * 3)))
                            for f = 1:3
                                if perturbed(ind_vinfi(f + (vind - 1) * 3)) < 0
                                    perturbed(ind_vinfi(f + (vind - 1) * 3)) = perturbed(ind_vinfi(f + (vind - 1) * 3)) + s_v * rand;
                                else
                                    perturbed(ind_vinfi(f + (vind - 1) * 3)) = perturbed(ind_vinfi(f + (vind - 1) * 3)) - s_v * rand;
                                end
                            end
                        end
                    end
                end
            elseif perturbed(j) > ub(j)
                perturbed(j) = ub(j);
                % Adjust vinfi if it equals vinff
                if N_flybys > 0
                    for vind = 1:N_flybys
                        if norm(perturbed(ind_vinfi((1:3) + (vind - 1) * 3))) == norm(perturbed(ind_vinff((1:3) + (vind - 1) * 3)))
                            for f = 1:3
                                if perturbed(ind_vinfi(f + (vind - 1) * 3)) < 0
                                    perturbed(ind_vinfi(f + (vind - 1) * 3)) = perturbed(ind_vinfi(f + (vind - 1) * 3)) + s_v * rand;
                                else
                                    perturbed(ind_vinfi(f + (vind - 1) * 3)) = perturbed(ind_vinfi(f + (vind - 1) * 3)) - s_v * rand;
                                end
                            end
                        end
                    end
                end
            end
        end

        % Re-optimize first phase
        consts(7) = 1;
        lbCurrent = lb(1:phaseSizes(1));
        ubCurrent = ub(1:phaseSizes(1));
        x = perturbed(1:phaseSizes(1));
        if whichThrust(1) == 1
            [optimised, ~, ~, output] = fmincon(@(x)obj_lofiSF(x, objInd(1)), x, A, b, Aeq, beq, lbCurrent, ubCurrent, @(x)con_lofiSF(x, consts), options);
        else
            [optimised, ~, ~, output] = fmincon(@(x)obj_lofiSF_coast(x, objInd(1)), x, A, b, Aeq, beq, lbCurrent, ubCurrent, @(x)con_lofiSF(x, consts), options);
        end
        perturbed(1:phaseSizes(1)) = optimised;

        % Re-optimize intermediate phases
        NpCurrent = 2;
        if Np > 2
            for i = 2:Np-1
                consts(7) = NpCurrent;
                x = perturbed(1:phaseSizes(i));
                lbCurrent = lb(1:phaseSizes(i));
                ubCurrent = ub(1:phaseSizes(i));
                if any(i == whichThrust)
                    [optimised, ~, ~, output] = fmincon(@(x)obj_lofiSF(x, objInd(i)), x, A, b, Aeq, beq, lbCurrent, ubCurrent, @(x)con_lofiSF(x, consts), options);
                else
                    if any(i <= whichThrust)
                        [optimised, ~, ~, output] = fmincon(@(x)obj_lofiSF(x, objInd(i)), x, A, b, Aeq, beq, lbCurrent, ubCurrent, @(x)con_lofiSF(x, consts), options);
                    else
                        [optimised, ~, ~, output] = fmincon(@(x)obj_lofiSF_coast(x, objInd(i)), x, A, b, Aeq, beq, lbCurrent, ubCurrent, @(x)con_lofiSF(x, consts), options);
                    end
                end
                perturbed(1:phaseSizes(i)) = optimised;
                NpCurrent = NpCurrent + 1;
            end
        end

        % Re-optimize final phase
        consts(7) = Np;
        x = perturbed;
        if whichThrust(end) == Np
            [optimised, m_optim, ~, output] = fmincon(@(x)obj_lofiSF(x, objInd(Np)), x, A, b, Aeq, beq, lb, ub, @(x)con_lofiSF(x, consts), options);
        else
            if any(whichThrust)
                [optimised, m_optim, ~, output] = fmincon(@(x)obj_lofiSF(x, objInd(Np)), x, A, b, Aeq, beq, lb, ub, @(x)con_lofiSF(x, consts), options);
            else
                [optimised, m_optim, ~, output] = fmincon(@(x)obj_lofiSF_coast(x, objInd(Np)), x, A, b, Aeq, beq, lb, ub, @(x)con_lofiSF(x, consts), options);
            end
        end

        % Save results to archive
        optim_archive(k, :) = optimised;
        m_archive(k) = m_optim;
        violation_archive(k) = output.constrviolation;
    end