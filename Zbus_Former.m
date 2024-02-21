clear all; 
clc;

[Z_bus, Y_bus] = Zbus_former('txt file path');


function [Z_bus, Y_bus] = Zbus_former(path)
    % Reading the data
    fid = fopen(path);
    %% BUS Data
    Line_String_Complete = fgetl(fid); 
    Bus_Data=[];
    Shunt_B = [];
    No_of_Buses=0;
    No_of_Shunt=0;
    while ischar(Line_String_Complete)
        Line_String_Complete=fgetl(fid);
        if(strcmp(Line_String_Complete(1:4),'-999')==1); 
            break; 
        end
        Line_String_Numeric=Line_String_Complete(15:end);
        Shunt_String_Numeric=Line_String_Complete(107:end);
        Line_Numeric=str2num(Line_String_Numeric);
        Shunt_Numeric=str2num(Shunt_String_Numeric);
        No_of_Buses=No_of_Buses+1;
        Bus_Data=[Bus_Data; [No_of_Buses Line_Numeric] ];
        Shunt_B = [Shunt_B; Shunt_Numeric] ;   
    end
    %% Line Data
    Line_String_Complete=fgetl(fid);
    Line_Data=[];
    No_of_Lines=0;
    while ischar(Line_String_Complete)
        Line_String_Complete=fgetl(fid);
    %     disp(['#' Line_String_Complete '#'] );
        if(strcmp(Line_String_Complete(1:4),'-999')==1);
            break;
        end
        Line_String_Numeric=Line_String_Complete(1:end);    
        Line_Numeric=str2num(Line_String_Numeric);
        No_of_Lines=No_of_Lines+1;
        Line_Data=[Line_Data; [No_of_Lines Line_Numeric]];
    end
    fclose(fid);
    % Data has been read
    
    % Processing the bus, shunt and line data
    from = Line_Data(:,2);
    to = Line_Data(:,3);
    to = round(to);
    R = Line_Data(:,8);
    X = Line_Data(:,9);
    B = Line_Data(:,10)*i;
    Line_Shunt_Z = imag(B.^-1)*i;  %to get rid of Inf values
    Z_line = R + X*i;
    number_of_busses = max(max(from), max(to));
    number_of_operations = length(from);
    % operation_order = Line_Data(:,1);
    Z_bus_not_redacted = zeros(number_of_busses+1,number_of_busses+1);
    old_buses = zeros(number_of_busses,1);
    shunts_from_lines = zeros(number_of_busses, 1);
    for line = 1:length(Line_Data)
        shunts_from_lines(Line_Data(line,2)) = Line_Data(line,10)/2;
        shunts_from_lines(Line_Data(line,3)) = Line_Data(line,10)/2;
    end
    external_Shunts = Shunt_B(:,2)*i;
    shunts_from_lines = shunts_from_lines*i;
    Shunt_Y = external_Shunts+shunts_from_lines;
    Shunt_Z = imag(Shunt_Y.^-1)*i;   %to get rid of Inf values
    % new_buses = zeros(number_of_busses,1);
    % Type one modification
    for bu = 1:length(Shunt_Z)
        if (Shunt_Z(bu)~=0 )
            old_buses(bu) = bu;
            Z_bus_not_redacted = type_one_modification(Z_bus_not_redacted, Shunt_Z, bu);
        end   
    end
    % Type two modification
    for b = 1:number_of_operations
        if xor((ismember(from(b), old_buses)==1), (ismember(to(b), old_buses)==1))
            if (ismember(from(b), old_buses)==1)
                old_one = from(b);
                new_one = to(b);
            else
                old_one = to(b);
                new_one = from(b);
            end
            old_buses(old_one) = old_one;
            old_buses(new_one) = new_one;
            Z_bus_not_redacted = type_two_modification(Z_bus_not_redacted, new_one, old_one, b, Z_line);
        end
    end
    % Type three modification
    for c = 1:number_of_operations %Reference bus is decided as bus 0
        if ((from(c)==0) | (to(c)==0)) & ((ismember(from(c), old_buses)==1) | (ismember(to(c), old_buses)==1))
            if (ismember(from(c), old_buses)==1) & (from(c)~=0)
                Z_bus_not_redacted = type_three_modification(Z_bus_not_redacted, from(c), c, Z_line, number_of_busses)
            end
            if (ismember(to(c), old_buses)==1) & (to(c)~=0)
                Z_bus_not_redacted = type_three_modification(Z_bus_not_redacted, to(c), c, Z_line, number_of_busses)
            end
        end
    end
    % Type four modification
    for d = 1:number_of_operations
        if ((ismember(from(d), old_buses)==1)) & ((ismember(to(d), old_buses)==1))
            Z_bus_not_redacted = type_four_modification(Z_bus_not_redacted, from(d), to(d), d, Z_line, number_of_busses);
        end
    end
    % Last edits to Z bus and forming the Y_bus
    Z_bus = Z_bus_not_redacted(1:number_of_busses, 1:number_of_busses);
    Y_bus_not_redacted = inv(Z_bus);
    Y_bool = abs(Y_bus_not_redacted)<1e-9;
    Y_bus_not_redacted_copy = Y_bus_not_redacted;
    Y_bus_not_redacted_copy(Y_bool) = 0;
    Y_bus = Y_bus_not_redacted_copy;
    
    % Functions for modifications
    function bus_matrix = type_one_modification(bus_matrix, shunt_data_matrix, shunt_bus_no);
            bus_matrix(shunt_bus_no, shunt_bus_no) = shunt_data_matrix(shunt_bus_no);
    end
    
    function bus_matrix = type_two_modification(bus_matrix, new, old, operation, impedances);
        for a = 1:new
            bus_matrix(new,a) = bus_matrix(old,a);
            bus_matrix(a, new) = bus_matrix(a,old);
            bus_matrix(new, new) = bus_matrix(old,old) + impedances(operation);
        end
    end
    
    function bus_matrix = type_three_modification(bus_matrix, old, operation, impedances, number_of_buses)
        payda = bus_matrix(old, old) + impedances(operation);
        for j = 1: number_of_buses
            for k = 1:number_of_buses
                bus_matrix(j, k) = bus_matrix(j, k) - (bus_matrix(old, k)*bus_matrix(j, old)/payda);
            end
        end
    end
    
    function bus_matrix = type_four_modification(bus_matrix, old_1, old_2, operation, impedances, number_of_buss)
        paydas = impedances(operation) + bus_matrix(old_1, old_1) + bus_matrix(old_2, old_2) - 2*bus_matrix(old_2, old_1);
        bus_matrix(:, number_of_buss+1) = bus_matrix(:,old_1)-bus_matrix(:,old_2);
        bus_matrix(number_of_buss+1, :) = bus_matrix(old_1,:)-bus_matrix(old_2,:);
        bus_matrix(number_of_buss+1, number_of_buss+1) = paydas;
        for j = 1: number_of_buss
            for k = 1:number_of_buss
                bus_matrix(j, k) = bus_matrix(j, k) - (bus_matrix(number_of_buss+1, k)*bus_matrix(j, number_of_buss+1)/paydas);
            end
        end
    end
end