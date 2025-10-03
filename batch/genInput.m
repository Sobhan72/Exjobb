%% Batch
E = 210e9;
v = 0.4;

sigy_base = [360e6, 300e6, 250e6];
E_base    = [E, 0.7*E, 0.5*E];
v_base    = [v, 0.5*v, 0.5*v];
sigc_values = [0.01, 0.1, 0.3, 0.5];

% x_values   = [0, 0, 0.7, 0.8, 0.9, 1.0];
x_values   = [0, 0.8, 1.0];
disp_values = [-1e-3, -9e-4, -8e-4];

n_cases = 3;
count = 0;

%ANISOTROPIC CASES
for c = 1:length(sigc_values)
    params.sigc = sigc_values(c);
for i = 1:length(x_values)
    x = x_values(i);

    for case_num = 1:n_cases
        % Rotate parameters according to case
        idx = mod((0:2) + case_num - 1, 3) + 1;

        for d = 1:length(disp_values)
            disp_val = disp_values(d);

            params.E1 = E_base(idx(1));
            params.E2 = E_base(idx(2));
            params.E3 = E_base(idx(3));

            if case_num == 1
                params.v12 = v_base(1);
                params.v13 = v_base(2);
                params.v23 = v_base(3);
            elseif case_num == 2
                params.v12 = v_base(2)*(E_base(3) / E_base(1)); % v31
                params.v13 = v_base(3)*(E_base(3) / E_base(2)); % v32
                params.v23 = v_base(1);
            else
                params.v12 = v_base(3);
                params.v13 = v_base(1)*(E_base(2) / E_base(1)); % v21
                params.v23 = v_base(2)*(E_base(3) / E_base(1)); % v31
            end

            params.sigy01 = sigy_base(idx(1));
            params.sigy02 = sigy_base(idx(2));
            params.sigy03 = sigy_base(idx(3));

            params.disp = disp_val;


            params.saveName = sprintf("Anisotrop_x=%02d_case%d_disp=%.0fe-4", ...
                round(x*10), case_num, -disp_val*1e4);

            filename = sprintf("input%d.mat", count);
            save(filename, 'x', 'params');

            count = count + 1;


        end
    end
end
end

%ISOTROPIC CASES
clearvars -except count E x_values disp_values sigc_values

v = 0.3;
params.sigy01 = 360e6;
params.sigy02 = 360e6;
params.sigy03 = 360e6;
params.E1 = E;
params.E2 = E;
params.E3 = E;
params.v12 = v;
params.v13 = v;
params.v23 = v;

for c = 1:length(sigc_values)
    params.sigc = sigc_values(c);
for i = 1:length(x_values)
    x = x_values(i);

    for d = 1:length(disp_values)
        disp_val = disp_values(d);
        params.disp = disp_val;

        params.saveName = sprintf("Isotrop_x=%02d_disp=%.0fe-4", ...
            round(x*10), -disp_val*1e4);

        filename = sprintf("input%d.mat", count);
        save(filename, 'x', 'params');

        count = count + 1;

    end
end
end


% %% Batch 
% E = 210e9;
% v = 0.3;
% % sigy_base = [360e6, 320e6, 280e6]; %01
% % E_base = [E, 0.85*E, 0.7*E];       %01
% % v_base = [v, 0.5*v, 0.5*v];        %01
% sigy_base = [360e6, 300e6, 250e6]; 
% E_base = [E, 0.7*E, 0.5*E];        
% v_base = [v, 0.5*v, 0.5*v];        
% 
% x_values = [0, 0, 0.7, 0.8, 0.9, 1.0];
% disp_values = [-1e-3, -9e-4, -8e-4];
% 
% n_cases = 3;
% 
% count = 0;
% 
% for i = 1:length(x_values)
%     x = x_values(i);    
% 
%     for case_num = 1:n_cases
%         % Rotate parameters according to case
%         idx = mod((0:2) + case_num - 1, 3) + 1;
% 
%         for d = 1:length(disp_values)
%             disp_val = disp_values(d);
% 
%             params.E1 = E_base(idx(1));
%             params.E2 = E_base(idx(2));
%             params.E3 = E_base(idx(3));
% 
%         if case_num == 1
%             params.v12 = v_base(1);
%             params.v13 = v_base(2);
%             params.v23 = v_base(3);
%         elseif case_num == 2
%             params.v12 = v_base(2)*(E_base(3) / E_base(1)); %v31 
%             params.v13 = v_base(3)*(E_base(3) / E_base(2)); %v32   
%             params.v23 = v_base(1);
%         else
%             params.v12 = v_base(3);    
%             params.v13 = v_base(1)*(E_base(2) / E_base(1)); %v21
%             params.v23 = v_base(2)*(E_base(3) / E_base(1)); %v31
%         end
% 
%             params.sigy01 = sigy_base(idx(1));
%             params.sigy02 = sigy_base(idx(2));
%             params.sigy03 = sigy_base(idx(3));
% 
%             params.disp = disp_val;
%             params.saveName = sprintf("Anisotrop_x=%02d_case%d_disp=%.0fe-4", ...
%                 round(x*10), case_num, -disp_val*1e4);
% 
%             filename = sprintf("input%d.mat", count);
%             save(filename, 'x', 'params');
% 
%             count = count + 1;
%         end
%     end
% end
% 
% clearvars -except count
% E = 210e9;
% v = 0.3;
% sigy_base = [360e6, 360e6,360e6]; 
% params.E1 = E;
% params.E2 = E;
% params.E3 = E;
% params.v12 = v;
% params.v13 = v;
% params.v23 = v;
% 
% x_values = [0, 0, 0.7, 0.8, 0.9, 1.0];
% disp_values = [-1e-3, -9e-4, -8e-4];
% 
% for i = 1:length(x_values)
%     x = x_values(i);    
% 
%          for d = 1:length(disp_values)
%             disp_val = disp_values(d);
%             params.disp = disp_val;
%             params.saveName = sprintf("Isotrop_x=%02d_disp=%.0fe-4", ...
%                 round(x*10), -disp_val*1e4);
% 
%             filename = sprintf("input%d.mat", count);
%             save(filename, 'x', 'params');
% 
%             count = count + 1;
%         end
% end
   

    % E = 210e9;
% v = 0.3;
% sigy_base = [360e6, 300e6, 250e6];
% E_base = [E, 0.7*E, 0.5*E];
% v_base = [v, 0.5*v, 0.5*v];
% 
% x_values = [0.7, 0.8, 0.9, 1.0];
% n_cases = 3;
% 
% count = 0;
% 
% fprintf("Jag är i %", pwd);
% fprintf("Börjar loopa!");
% for i = 1:length(x_values)
%     x = x_values(i);    
% 
%     for case_num = 1:n_cases
%         % Rotate parameters according to case
%         idx = mod((0:2) + case_num - 1, 3) + 1;
% 
%         params.E1 = E_base(idx(1));
%         params.E2 = E_base(idx(2));
%         params.E3 = E_base(idx(3));
% 
%         if case_num == 1
%             params.v12 = v_base(1);
%             params.v13 = v_base(2);
%             params.v23 = v_base(3);
%         elseif case_num == 2
%             params.v12 = v_base(2)*(E_base(3) / E_base(1)); %v31 
%             params.v13 = v_base(3)*(E_base(3) / E_base(2)); %v32   
%             params.v23 = v_base(1);
%         else
%             params.v12 = v_base(3);    
%             params.v13 = v_base(1)*(E_base(2) / E_base(1)); %v21
%             params.v23 = v_base(2)*(E_base(3) / E_base(1)); %v31
%         end
% 
%         params.sigy01 = sigy_base(idx(1));
%         params.sigy02 = sigy_base(idx(2));
%         params.sigy03 = sigy_base(idx(3));
% 
%         % Save name and filename
%         params.saveName = sprintf("OptDesign_batch5_x=%02d_case%d", round(x*10), case_num);
%         filename = sprintf("input%d.mat", count);
%         save(filename, 'x', 'params');
% 
%         fprintf("Sparat inputs till fil %s", filename);
%         count = count + 1;
%     end
% end
% fprintf("Färdig!")



% %% Batch3
% File: generateInputs.m

% fileIdx = 0;
% 
% for xi = 1:10                    % indexera x
%     x = 0.1 * xi;                
%     for pq = 0:1                 % params.PQ can be 0 or 1
%         params = struct();       % start fresh for each file
%         params.PQ = pq;          
%         if pq == 1
%             params.p = 1.5;
%             params.q = 1;
%         else
%             params.p = 3;
%             params.q = 2;
%         end
% 
%         params.saveName = sprintf("OptDesign_batch3_x=%02d_PQ=%d", round(x*10), pq);
% 
%         filename = sprintf('input%d.mat', fileIdx);
% 
%         save(filename, 'x', 'params');
% 
%         fileIdx = fileIdx + 1; 
%     end
% end




% %% Batch 2
% for caseNum = 1:4
%     switch caseNum
%         case 1
%             params.rampB = 1;
%             params.rampPQ = false;
%             params.p = 3;
%             params.q = 2.5;
%             params.beta = 1;
%         case 2
%             params.rampB = 2;
%             params.rampPQ = false;
%             params.p = 3;
%             params.q = 2.5;
%             params.beta = 1e-6;
%         case 3
%             params.rampB = 1;
%             params.rampPQ = true;
%             params.p = 1.5;
%             params.q = 1;
%             params.beta = 1;
%         case 4
%             params.rampB = 2;
%             params.rampPQ = true;
%             params.p = 1.5;
%             params.q = 1;
%             params.beta = 1e-6;
%     end
% 
%     re_values = [3, 2];
%     disp_values = [-1.5e-3, -1e-3, -5e-4];
% 
%     index = (caseNum - 1) * 6;
% 
%     for i = 1:length(re_values)
%         for j = 1:length(disp_values)
%             params.re = re_values(i);
%             params.disp = disp_values(j);
% 
%             params.saveName = sprintf("OptDesign_batch2_B=%d_PQ=%d_re=%d_disp=%g", ...
%                                       params.rampB, params.rampPQ, params.re, params.disp);
% 
%             x = [];
%             filename = sprintf('input%d.mat', index);
%             save(filename, 'x', 'params');
%             index = index + 1;
%         end
%     end
% end