function [optimised, m_optim, constrviolation] = basinhop(k, optimised, sigmas, MBH_tail, ...
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


        % Ensure all variables are local to avoid issues in parfor
        probSize = length(optimised);
        constsCopy = consts; % Avoid broadcast issues. TERRIBLE
        m_optim = 0; % Initialize m_optim
        output = struct(); % Initialize output as an empty struct

        % Perturb all the decision variables
        pm = fix(rand(1, probSize) + 0.5);
        pm(~pm) = -1;

        % TODO: This doesn't need to be done twice.
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
        % TODO: Must be better way to do this!
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
        % TODO: See no reason this cant be in the loop as index 1
        constsCopy(7) = 1;
        lbCurrent = lb(1:phaseSizes(1));
        ubCurrent = ub(1:phaseSizes(1));
        x = perturbed(1:phaseSizes(1));
        if whichThrust(1) == 1
            [optimised, ~, ~, ~] = fmincon(@(x)obj_lofiSF(x, objInd(1)), x, A, b, Aeq, beq, lbCurrent, ubCurrent, @(x)con_lofiSF(x, constsCopy), options);
        else
            [optimised, ~, ~, ~] = fmincon(@(x)obj_lofiSF_coast(x, objInd(1)), x, A, b, Aeq, beq, lbCurrent, ubCurrent, @(x)con_lofiSF(x, constsCopy), options);
        end
        perturbed(1:phaseSizes(1)) = optimised;

        % Re-optimize intermediate phases
        NpCurrent = 2;
        for i = 2:Np-1
            constsCopy(7) = NpCurrent;
            x = perturbed(1:phaseSizes(i));
            lbCurrent = lb(1:phaseSizes(i));
            ubCurrent = ub(1:phaseSizes(i));

            % TODO: Simplify this
            if any(i == whichThrust)
                [optimised, ~, ~, ~] = fmincon(@(x)obj_lofiSF(x, objInd(i)), x, A, b, Aeq, beq, lbCurrent, ubCurrent, @(x)con_lofiSF(x, constsCopy), options);
            else
                if any(i <= whichThrust)
                    [optimised, ~, ~, ~] = fmincon(@(x)obj_lofiSF(x, objInd(i)), x, A, b, Aeq, beq, lbCurrent, ubCurrent, @(x)con_lofiSF(x, constsCopy), options);
                else
                    [optimised, ~, ~, ~] = fmincon(@(x)obj_lofiSF_coast(x, objInd(i)), x, A, b, Aeq, beq, lbCurrent, ubCurrent, @(x)con_lofiSF(x, constsCopy), options);
                end
            end
            perturbed(1:phaseSizes(i)) = optimised;
            NpCurrent = NpCurrent + 1;
        end

        % Re-optimize final phase
        constsCopy(7) = Np;
        x = perturbed;
        % TODO: 
        % surely `whichThrust(end) == Np` is a subset of `any(whichThrust)`
        % This can also be put into loop.
        if whichThrust(end) == Np
            [optimised, m_optim, ~, output] = fmincon(@(x)obj_lofiSF(x, objInd(Np)), x, A, b, Aeq, beq, lb, ub, @(x)con_lofiSF(x, constsCopy), options);
        else
            if any(whichThrust)
                [optimised, m_optim, ~, output] = fmincon(@(x)obj_lofiSF(x, objInd(Np)), x, A, b, Aeq, beq, lb, ub, @(x)con_lofiSF(x, constsCopy), options);
            else
                [optimised, m_optim, ~, output] = fmincon(@(x)obj_lofiSF_coast(x, objInd(Np)), x, A, b, Aeq, beq, lb, ub, @(x)con_lofiSF(x, constsCopy), options);
            end
        end
        
        constrviolation = output.constrviolation;
        % Where is optimised compared with

        % TODO: 'optimised', 'm_optim' and 'output' value being returned is just for Final Phase? Is
        % this correct?
    end