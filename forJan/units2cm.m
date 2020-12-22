function position_cm = units2cm(position_units)
%UNITS2CM converts Jan Sigurd's unit data to cm

position_cm = [position_units(:,1), position_units(:,2:5).*100;];


end

