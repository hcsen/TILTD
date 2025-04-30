function [optimised, m_optim, constrviolation] = basinhop(k, best, sigmas, MBH_tail, ...
        MBH_theta, rho_hop, t0Hop, dtHop, dtIndex, lb, ub, ind_vinfi, ind_vinff, s_v, N_flybys, ...
        phaseSizes, objInd, whichThrust, consts, options, A, b, Aeq, beq, Np)


        % k = index
        % optimised = reduction variable
        % sigmas = broadcast variable
        % MBH_tail = broadcast variable
        % MBH_theta = broadcast variable
        % rho_hop = broadcast variable
        % t0Hop = broadcast variable

        % dtHop = broadcast variable
        % dtIndex = broadcast variable
        % lb = broadcast variable
        % ub = broadcast variable
        % ind_vinfi = broadcast variable
        % ind_vinff = broadcast variable
        % s_v = broadcast variable
        % N_flybys = broadcast variable
        % phaseSizes = broadcast variable
        % objInd = broadcast variable
        % whichThrust = broadcast variable
        % constsCopy = INVALID
        % options = broadcast variable
        % A = broadcast variable
        % b = broadcast variable
        % Aeq = broadcast variable
        % beq = broadcast variable
        % Np = broadcast variable

        % Reseed random based on index;
        baserng = rng();
        rng(baserng.Seed + k);

        % Ensure all variables are local to avoid issues in parfor
        constsCopy = consts; % Avoid broadcast issues. TERRIBLE
        m_optim = 0; % Initialize m_optim
        output = struct(); % Initialize output as an empty struct
        

        probSize = length(best);

        % Perturb all the decision variables
        pm = fix(rand(1, probSize) + 0.5);
        pm(~pm) = -1;

        % TODO: This doesn't need to be done twice.
        perturbed = best + gprnd(MBH_tail, sigmas, MBH_theta * ones(1, probSize)) .* pm;

        % Hop to a new random guess if needed
        if rand < rho_hop
            disp('Random hop');
            perturbed(1) = perturbed(1) + 2 * t0Hop * rand - t0Hop; % Hop the launch epoch
            for j = 1:Np
                perturbed(dtIndex(j)) = perturbed(dtIndex(j)) + 2 * dtHop * rand - dtHop; % Hop each phase tof
            end
        end

        % Ensure decision variables are within bounds
        % TODO: Must be better way to do this!
        for j = 1:probSize
            if perturbed(j) < lb(j)
                perturbed(j) = lb(j);
                % Adjust vinfi if it equals vinff
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
            elseif perturbed(j) > ub(j)
                perturbed(j) = ub(j);
                % Adjust vinfi if it equals vinff
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

        % Re-optimize phases
        for i = 1:Np
            fprintf("MBH Hop (%u).... Phase(%u/%u)\n",k, i, Np);

            constsCopy(7) = i;
            x = perturbed(1:phaseSizes(i));
            lbCurrent = lb(1:phaseSizes(i));
            ubCurrent = ub(1:phaseSizes(i));

            if any(i == whichThrust)
                [optimised, ~, ~, output] = fmincon(@(x)obj_lofiSF(x, objInd(i)), x, A, b, Aeq, beq, lbCurrent, ubCurrent, @(x)con_lofiSF(x, constsCopy), options);
            else
                [optimised, ~, ~, output] = fmincon(@(x)obj_lofiSF_coast(x, objInd(i)), x, A, b, Aeq, beq, lbCurrent, ubCurrent, @(x)con_lofiSF(x, constsCopy), options);
            end
            perturbed(1:phaseSizes(i)) = optimised;
        end
        fprintf("MBH Hop (%u).... Done\n",k);
        constrviolation = output.constrviolation;
    end