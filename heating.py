def thermal_analysis(t, dt, w, nx, u, therm_con, spec_heat, density):
    # t - time(s)
    # dt - timestep
    # w - width (m)
    # nx - number of spatial steps

    nt = int(t / dt)  # nt - number of time steps
    dx = w / nx  # dx - spatial step

    # Extend data for small timesteps
    x_given = np.arange(len(u))
    x_interpolated = np.linspace(0, len(u) - 1, nt)
    y_interpolated = np.interp(x_interpolated, x_given, u)
    interpolated_data = np.column_stack((x_interpolated, y_interpolated))

    # Calculate alpha and p from material properties
    alpha = therm_con / (spec_heat * density)
    p = alpha * dt / dx ** 2
    data = np.zeros((nt, nx))

    # Append initial condition to data array
    data[:, 0] = interpolated_data[:, 1]
    data[0, :] = u[0]
    time_array = np.linspace(0, t, nt)

    # Forward differencing calculations
    for x in range(nt - 1):
        x = x + 1
        for n in range(nx - 1):
            n = n + 1
            if n == nx - 1:  # For boundary condition
                data[x, n] = (1 - 2 * p) * data[x - 1, n] + p * (2 * data[x - 1, n - 1])  # Neumann boundary
            else:
                data[x, n] = (1 - 2 * p) * data[x - 1, n] + p * (data[x - 1, n - 1] + data[x - 1, n + 1])

    np.set_printoptions(formatter={'float': '{:.4g}'.format})

    leth = []
    for i in range(nx):
        l_ = lethality(3, 10, dt, data[:, i])  # Calculate cycle lethality at each spatial step
        # l_ = lethality(3, 10, dt, data[:, i])      - lethality for standard IO
        # l_ = lethality(4.88, 11.1, dt, data[:, i]) - lethality for worst-case IO

        leth.append(l_)
    # Calculate the smallest lethality and spatial step
    smallest_lethality = min(leth)
    print("All lethalities: ")
    print(leth)
    depth = (np.argmin(leth)) * (w / (nx-1))
    print(f"The lowest lethality is: {smallest_lethality} at a depth of {depth}m")

    # Print temperature curve at each spatial steps
    for col in data.T:
        plt.plot(time_array, col)

    plt.xlabel('Time(s)')
    plt.ylabel('Temperature(Celsius)')
    plt.title('Temperature Distribution')
    plt.ylim([0, 150])  # Set y-axis limits to 0-150
    plt.grid(True)
    plt.show()

    return smallest_lethality


# Calculate lethality
def lethality(d_value, z_value, nt, temp):

    TDT_array = []
    for i in temp:
        l = d_value * 10 ** ((121 - i) / z_value)
        TDT_array.append(l)  # Append Thermal Death Times

    lethality_rate_array = np.power(TDT_array, -1)  # Calculate lethality rate
    lethality_array = lethality_rate_array * nt * 1 / 60
    lethality = np.sum(lethality_array)
    print("Total lethality")
    print(lethality)
    return lethality


# Temperature data for 121 cycle
Heat_121 = [10, 50, 80, 110, 115, 121]
Hold_121 = [119.5, 121.5]
Cool_121 = [120, 70, 50, 40, 30]

# Temperature data for 127 cycle
Heat_127 = [10, 50, 80, 110, 115, 127]
Hold_127 = [126.5, 127.5]
Cool_127 = [120, 70, 50, 40, 30]

# Temperature data for 134 cycle
Heat_134 = [10, 50, 80, 110, 115, 135]
Hold_134 = [133.5, 134.5]
Cool_134 = [120, 70, 50, 40, 30]


# Titanium - hold cycle

# leth_ti_121_hold = thermal_analysis(40.1*60, 0.01, 0.1, 5, Hold_121, 11.4, 0.5223, 4540) # Holding
# print(f"Total Lethality: {leth_ti_121_hold}")
# leth_ti_127_hold = thermal_analysis(25.3*60, 0.01, 0.1, 5, Hold_127, 11.4, 0.5223, 4540) # Holding
# print(f"Total Lethality: {leth_ti_127_hold}")
# leth_ti_134_hold = thermal_analysis(5.9 * 60, 0.03, 0.1, 5, Hold_134, 11.4, 0.5223, 4540)  # Holding
# print(f"Total Lethality: {leth_ti_134_hold}")


# Titanium - entire cycle

# leth_ti_121_heat = thermal_analysis(30*60, 0.01, 0.1, 5, Heat_121, 11.4, 0.5223, 4540) # Heating
# leth_ti_121_hold = thermal_analysis(58*60, 0.01, 0.1, 5, Hold_121, 11.4, 0.5223, 4540) # Holding
# leth_ti_121_cool = thermal_analysis(60*60, 0.01, 0.1, 5, Cool_121, 11.4, 0.5223, 4540) # Cooling
# leth_ti_121_total = leth_ti_121_heat + leth_ti_121_hold + leth_ti_121_cool
# print(leth_ti_121_total)
# print(f"Total Lethality: {leth_ti_121_total}")

# leth_ti_127_heat = thermal_analysis(30*60, 0.01, 0.1, 5, Heat_127, 11.4, 0.5223, 4540) # Heating
# leth_ti_127_hold = thermal_analysis(11.1*60, 0.01, 0.1, 5, Hold_127, 11.4, 0.5223, 4540) # Holding
# leth_ti_127_cool = thermal_analysis(60*60, 0.01, 0.1, 5, Cool_127, 11.4, 0.5223, 4540) # Cooling
# leth_ti_127_total = leth_ti_127_heat + leth_ti_127_hold + leth_ti_127_cool
# print(leth_ti_127_total)

# leth_ti_134_heat = thermal_analysis(10, 0.01, 0.03, 5, Cool_121, 11.4, 0.5223, 4540) # Heating
# leth_ti_134_hold = thermal_analysis(10*60, 0.03, 0.1, 5, Hold_134, 11.4, 0.5223, 4540) # Holding
# leth_ti_134_cool = thermal_analysis(60*60,0.03, 0.1, 5, Cool_134, 11.4, 0.5223, 4540) # Cooling
# leth_ti_134_total = leth_ti_134_heat + leth_ti_134_hold + leth_ti_134_cool
# print(leth_ti_134_total)


# Polycarbonate - entire cycle

# leth_ca_121_heat = thermal_analysis(30*60, 0.01, 0.03, 5, Heat_121, 0.23, 1.2, 1200) # Heating
# leth_ca_121_hold = thermal_analysis(90.5*60, 0.01, 0.03, 5, Hold_121, 0.23, 1.2, 1200) # Holding
# leth_ca_121_cool = thermal_analysis(60*60, 0.01, 0.03, 5, Cool_121, 0.23, 1.2, 1200) # Cooling
# leth_ca_121_total = leth_ca_121_heat + leth_ca_121_hold + leth_ca_121_cool
# print(leth_ca_121_total)
# print(f"Total Lethality: {leth_ca_121_total}")

# leth_ca_127_heat = thermal_analysis(30*60, 0.01, 0.1, 5, Heat_127, 0.23, 1.2, 1200) # Heating
# leth_ca_127_hold = thermal_analysis(11.1*60, 0.01, 0.1, 5, Hold_127, 0.23, 1.2, 1200) # Holding
# leth_ca_127_cool = thermal_analysis(60*60, 0.01, 0.1, 5, Cool_127, 0.23, 1.2, 1200) # Cooling
# leth_ca_127_total = leth_ti_127_heat + leth_ti_127_hold + leth_ti_127_cool
# print(leth_ti_127_total)

# leth_ca_134_heat = thermal_analysis(10, 0.01, 0.03, 5, Cool_121, 0.23, 1.2, 1200) # Heating
# leth_ca_134_hold = thermal_analysis(10*60, 0.03, 0.1, 5, Hold_134, 0.23, 1.2, 1200) # Holding
# leth_ca_134_cool = thermal_analysis(60*60,0.03, 0.1, 5, Cool_134, 0.23, 1.2, 1200) # Cooling
# leth_ca_134_total = leth_ca_134_heat + leth_ca_134_hold + leth_ca_134_cool
# print(leth_ca_134_total)

