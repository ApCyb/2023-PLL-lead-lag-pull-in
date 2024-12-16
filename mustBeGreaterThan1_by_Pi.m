function mustBeGreaterThan1_by_Pi(k)
% Custom validation: k > 1/pi
if k <= 1/pi
    error("k must be greater than 1/pi.");
end
end