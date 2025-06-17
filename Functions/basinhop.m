function [optimised, m_optim, constrviolation] = basinhop(k, best, seed, sigmas, MBH_tail, ...
            rho_hop, t0Hop, dtHop, dtIndex, lb, ub, phaseSizes, objInd, whichThrust, consts, ...
            options, A, b, Aeq, beq, Np)

        % Reseed from input param.
        rng(seed);

        % Ensure all variables are local to avoid issues in parfor
        constsCopy = consts; % Avoid broadcast issues. TERRIBLE
        m_optim = 0; % Initialize m_optim
        output = struct(); % Initialize output as an empty struct
        

        probSize = length(best);

        % Perturb all the decision variables
        
        % New method using std
        % range = ub - lb;
        % Change 0.25 (2 standard deviations) to variable.
        % Keeping hardcoded now to avoid changing config.
        % perturbed = normrnd(best, 0.25 * range); 

        
        % Old method.
        pm = fix(rand(1, probSize) + 0.5);
        pm(~pm) = -1;
        perturbed = best + gprnd(MBH_tail, sigmas, 0) .* pm;
        fprintf("Worker Seed: %u\n", rng().Seed);
        irand = rand;

        % Hop to a new random guess if needed
        if  irand < rho_hop
            fprintf('Hop (%f < %f)\n', irand, rho_hop);
            % What is this doing.
            perturbed(1) = perturbed(1) + 2 * t0Hop * rand - t0Hop; % Hop the launch epoch
            for j = 1:Np
                perturbed(dtIndex(j)) = perturbed(dtIndex(j)) + 2 * dtHop * rand - dtHop; % Hop each phase tof
            end
        else
            fprintf('No Hop (%f > %f)\n', irand, rho_hop);
        end

        outsideboundcount = 0;

        % Ensure decision variables are within bounds
        % TODO: Must be better way to do this!
        for j = 1:probSize
            if perturbed(j) < lb(j)
                %fprintf('Perturbed Value Outside Bounds! %.4e < %.4e\n', perturbed(j),  lb(j));
                outsideboundcount = outsideboundcount + 1;
                perturbed(j) = lb(j);
            elseif perturbed(j) > ub(j)
                %fprintf('Perturbed Value Outside Bounds! %.4e > %.4e\n', perturbed(j),  ub(j));
                outsideboundcount = outsideboundcount + 1;
                perturbed(j) = ub(j);
            end
        end
        
        if outsideboundcount > 0
            fprintf('(%u / %u) perturbed values fell outside of bounds and were corrected.\n', outsideboundcount,  probSize);
        end

        % Re-optimize phases
        for i = 1:Np
            fprintf("MBH Step (%u).... Phase(%u / %u)\n",k, i, Np);

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
        fprintf("MBH Step (%u).... Done\n",k);
        constrviolation = output.constrviolation;
    end