function grid = explicit_uncertainty_grid(cfg)
%EXPLICIT_UNCERTAINTY_GRID Return the exact 15 Table-7 cases.
    if nargin < 1 || isempty(cfg)
        cfg = paperstudies.saito.config();
    end
    published = cfg.tables.table7.values;
    grid.columns = ["case_id", "range_to_go_km", "density_scale", "cd_scale"];
    grid.values = [published(:,1:2), published(:,3:4) / 100];
    grid.count = size(grid.values, 1);
    grid.source = "TABLE_7_PUBLISHED_EXACT";
    grid.random_sampling = false;
end
