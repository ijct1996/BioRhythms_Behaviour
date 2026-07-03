function [timeMinutes, TsMinutes] = hsub_infer_time(timeCol, timeName)
%HSUB_INFER_TIME Infer minute time base and sampling interval.
%
%   Matches legacy infer_time_minutes in harmonic_subtraction_v6.m.
%   Pipeline_Input uses Time (hr) with 10 min steps over ~5–10 d.

    if isdatetime(timeCol)
        t = timeCol(:);
        dt = minutes(diff(t));
        TsMinutes = median(dt, 'omitnan');
        timeMinutes = minutes(t - t(1));
        return;
    end

    if ~isnumeric(timeCol)
        error('hsub_infer_time:UnsupportedType', 'Time column must be numeric or datetime.');
    end

    t = double(timeCol(:));
    nm = lower(string(timeName));

    if contains(nm, "min")
        timeMinutes = t;
    elseif contains(nm, "hr") || contains(nm, "hour")
        timeMinutes = t * 60;
    elseif contains(nm, "day")
        timeMinutes = t * 24 * 60;
    else
        dtRaw = diff(t);
        medStep = median(dtRaw, 'omitnan');
        % Heuristic: long span + small steps → values are hours
        if max(t, [], 'omitnan') > 24 && medStep > 0 && medStep < 1
            timeMinutes = t * 60;
        else
            timeMinutes = t;
        end
    end

    dt = diff(timeMinutes);
    TsMinutes = median(dt, 'omitnan');
    if ~isfinite(TsMinutes) || TsMinutes <= 0
        error('hsub_infer_time:BadSampling', ...
            'Sampling interval could not be inferred from Time column.');
    end
end
