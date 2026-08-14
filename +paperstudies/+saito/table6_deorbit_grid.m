function grid = table6_deorbit_grid(cfg)
%TABLE6_DEORBIT_GRID Exact ignition epochs and published aggregate bounds.
%   The paper does not print the 27 individual entry states or their exact
%   range-to-go values. Those fields remain NaN instead of being fabricated
%   by linear interpolation between the published bounds.

    if nargin < 1 || isempty(cfg)
        cfg = paperstudies.saito.config();
    end
    source = cfg.tables.table6;
    offsets = source.deorbit_ignition_offset_s(:);
    count = numel(offsets);
    if count ~= source.total_number_of_cases
        error('paperstudies:saito:table6:CountMismatch', ...
              'Published Table-6 offsets must produce exactly 27 cases.');
    end

    start_epoch = datetime(char(source.deorbit_ignition_start_utc), ...
        'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z''', 'TimeZone', 'UTC');
    ignition_epochs = start_epoch + seconds(offsets);
    ignition_epochs.Format = 'yyyy-MM-dd''T''HH:mm:ss''Z''';

    grid.columns = ["case_id", "ignition_offset_s", ...
                    "individual_range_to_go_km"];
    grid.values = [(1:count).', offsets, nan(count,1)];
    grid.ignition_epoch_utc = string(ignition_epochs);
    grid.center_of_gravity_mm = source.center_of_gravity_mm;
    grid.elevation_angle_deg = source.elevation_angle_deg;
    grid.published_range_to_go_bounds_km = source.entry_range_to_go_bounds_km;
    grid.individual_range_to_go_available = false;
    grid.individual_entry_states_available = false;
    grid.count = count;
    grid.source = "TABLE_6_PUBLISHED_TIMES_AND_AGGREGATE_BOUNDS";
end
