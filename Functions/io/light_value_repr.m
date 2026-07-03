function v = light_value_repr(lightVec)
%LIGHT_VALUE_REPR Scalar light duration when constant; else NaN.
    try
        if isnumeric(lightVec) || islogical(lightVec)
            u = unique(lightVec(~isnan(lightVec)));
            if isscalar(u)
                v = u;
            else
                v = NaN;
            end
        else
            v = NaN;
        end
    catch
        v = NaN;
    end
end
