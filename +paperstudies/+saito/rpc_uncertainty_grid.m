function grid = rpc_uncertainty_grid(cfg)
%RPC_UNCERTAINTY_GRID Cartesian product of the exact Table-9 axes.
    if nargin < 1 || isempty(cfg)
        cfg = paperstudies.saito.config();
    end

    ranges = cfg.tables.table9.range_to_go_km;
    density = cfg.tables.table9.density_scale_percent / 100;
    cd = cfg.tables.table9.drag_coefficient_scale_percent / 100;
    guidance_names = cfg.tables.table9.guidance_approach;
    count = numel(ranges) * numel(density) * numel(cd) * numel(guidance_names);
    values = zeros(count, 5);
    guidance = strings(count, 1);

    row = 0;
    for range_index = 1:numel(ranges)
        for density_index = 1:numel(density)
            for cd_index = 1:numel(cd)
                for guidance_index = 1:numel(guidance_names)
                    row = row + 1;
                    values(row,:) = [row, ranges(range_index), ...
                        density(density_index), cd(cd_index), guidance_index];
                    guidance(row) = guidance_names(guidance_index);
                end
            end
        end
    end

    grid.columns = ["case_id", "range_to_go_km", "density_scale", ...
                    "cd_scale", "guidance_method_id"];
    grid.values = values;
    grid.guidance_method = guidance;
    grid.guidance_method_lookup = guidance_names;
    grid.count = count;
    grid.source = "TABLE_9_CARTESIAN_PRODUCT_PUBLISHED_EXACT";
    grid.random_sampling = false;
end
