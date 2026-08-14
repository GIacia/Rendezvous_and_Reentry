function ld_command = eq13_ld_command(ld_current, F1, R_current, R_guidance)
%EQ13_LD_COMMAND Saito Eq. (13), with caller-controlled units.
    ld_command = ld_current + F1 .* (R_current - R_guidance);
end
