function value = get_bool(s, name, default_value)
    value = get_field(s, name, default_value);
    if islogical(value) && isscalar(value)
        return;
    end
    if isnumeric(value) && isreal(value) && isscalar(value) && ...
            isfinite(value) && any(value == [0, 1])
        value = logical(value);
        return;
    end
    txt = lower(strtrim(string(value)));
    if isscalar(txt) && any(txt == ["true", "1", "yes", "on"])
        value = true;
    elseif isscalar(txt) && any(txt == ["false", "0", "no", "off"])
        value = false;
    else
        error('reentry_core:getBool:InvalidLogical', ...
            '%s must be a scalar logical value.', name);
    end
end
