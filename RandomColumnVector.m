function column_vector = RandomColumnVector(length, minValue, maxValue)
    column_vector = randi([minValue, maxValue], length, 1);
end