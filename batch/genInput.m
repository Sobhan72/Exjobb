%% Batch 1
for caseNum = 1:4
    switch caseNum
        case 1
            params.rampB = 1;
            params.rampPQ = false;
            params.p = 3;
            params.q = 2.5;
            params.beta = 1;
        case 2
            params.rampB = 2;
            params.rampPQ = false;
            params.p = 3;
            params.q = 2.5;
            params.beta = 1e-6;
        case 3
            params.rampB = 1;
            params.rampPQ = true;
            params.p = 1.5;
            params.q = 1;
            params.beta = 1;
        case 4
            params.rampB = 2;
            params.rampPQ = true;
            params.p = 1.5;
            params.q = 1;
            params.beta = 1e-6;
    end

    re_values = [3, 2];
    disp_values = [-1.5e-3, -1e-3, -5e-4];

    index = (caseNum - 1) * 6;

    for i = 1:length(re_values)
        for j = 1:length(disp_values)
            params.re = re_values(i);
            params.disp = disp_values(j);

            params.saveName = sprintf("OptDesign_batch2_B=%d_PQ=%d_re=%d_disp=%g", ...
                                      params.rampB, params.rampPQ, params.re, params.disp);

            x = [];
            filename = sprintf('input%d.mat', index);
            save(filename, 'x', 'params');
            index = index + 1;
        end
    end
end