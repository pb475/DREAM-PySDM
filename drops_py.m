function [output] = drops_py(input)
    persistent PySDM;
    persistent PySDM_environments;
    persistent PySDM_products;
    persistent PySDM_initialisation;
    persistent PySDM_spectra;
    persistent PySDM_spectral_sampling;
    persistent PySDM_physics;
    persistent PySDM_dynamics;
    persistent PySDM_backends;
    persistent si;
    if isempty(PySDM)
        PySDM = py.importlib.import_module('PySDM');
        PySDM_environments = py.importlib.import_module('PySDM.environments');
        PySDM_products = py.importlib.import_module('PySDM.products');
        PySDM_initialisation = py.importlib.import_module('PySDM.initialisation');
        PySDM_spectra = py.importlib.import_module('PySDM.initialisation.spectra');
        PySDM_spectral_sampling = py.importlib.import_module('PySDM.initialisation.sampling.spectral_sampling');
        PySDM_physics = py.importlib.import_module('PySDM.physics');
        PySDM_dynamics = py.importlib.import_module('PySDM.dynamics');
        PySDM_backends = py.importlib.import_module('PySDM.backends');
        si = PySDM_physics.constants.si;
    end
    formulae = PySDM.Formulae(pyargs( ...
        'constants', py.dict(pyargs( ...
            'sgm_w', input.sigma * si.joule / si.metre^2, ...  
            'MAC', input.MAC ...
        )) ...
    ));

    assert(length(input.kappa) == length(input.n_tot))
    assert(length(input.kappa) == length(input.meanr))
    assert(length(input.kappa) == length(input.gstdv))

    r_dry = ones(1, input.n_bins, 'double') * NaN;
    n_vol = ones(1, input.n_bins, 'double') * NaN;
    kappa = ones(1, input.n_bins, 'double') * NaN;

    counts = zeros(1, length(input.kappa), 'int32');
    for k=1:length(counts)
        counts(k) = floor(input.n_bins / length(counts));
    end

    if length(counts) > 1    
        offsets = circshift(cumsum(counts), 1);
        counts(length(counts)) = counts(length(counts)) + input.n_bins - offsets(1);
        offsets(1) = 0;
    else
        offsets = [0];
    end

    for k=1:length(input.kappa)
        n_tot = input.n_tot{k};
        meanr = input.meanr{k};
        gstdv = input.gstdv{k};
        
        assert(length(n_tot) == length(meanr))
        assert(length(gstdv) == length(meanr))
        if length(n_tot) == 1
            spectrum = PySDM_spectra.Lognormal(pyargs(...
                'norm_factor', n_tot(1), ...
                'm_mode', meanr(1), ...
                's_geom', gstdv(1) ...
            ));
        elseif length(n_tot) == 2
            spectrum = PySDM_spectra.Sum(py.list({
                PySDM_spectra.Lognormal(pyargs(...
                    'norm_factor', n_tot(1), ...
                    'm_mode', meanr(1), ...
                    's_geom', gstdv(1) ...
                )), ...
                PySDM_spectra.Lognormal(pyargs(...
                    'norm_factor', n_tot(2), ...
                    'm_mode', meanr(2), ...
                    's_geom', gstdv(2) ...
                )), ...        
            }));
        else
            assert(false)
        end
        tmp = PySDM_spectral_sampling.ConstantMultiplicity(pyargs(...
            'spectrum', spectrum...
        )).sample(pyargs('n_sd', counts(k)));
        r_dry(offsets(k)+1 : offsets(k) + counts(k)) = double(py.numpy.ascontiguousarray(tmp{1}));
        n_vol(offsets(k)+1 : offsets(k) + counts(k)) = double(tmp{2});
        kappa(offsets(k)+1 : offsets(k) + counts(k)) = input.kappa{k};
    end
    r_dry = py.numpy.array(r_dry);
    n_vol = py.numpy.array(n_vol);
    kappa = py.numpy.array(kappa);

    pv0 = input.RH * formulae.saturation_vapour_pressure.pvs_Celsius(input.T - PySDM_physics.constants.T0);
    q0 = PySDM_physics.constants_defaults.eps * pv0 / (input.p - pv0);
    R = (1 + q0) * (PySDM_physics.constants_defaults.Rv / (1 / q0 + 1) + PySDM_physics.constants_defaults.Rd / (1 + q0));
    rho0 = input.p / input.T / R;
    mass_of_air = 1000 * si.kg;
    environment = PySDM_environments.Parcel(pyargs( ...
        'dt', input.dt * si.s, ...
        'mass_of_dry_air', mass_of_air, ...
        'p0', input.p, ...
        'q0', q0, ...
        'T0', input.T * si.K, ...
        'w', input.w * si.m / si.s ...
    ));

    builder = PySDM.Builder(pyargs( ...
        'backend', PySDM_backends.CPU(formulae), ...
        'n_sd', int32(input.n_bins) ...
    ));
    builder.set_environment(environment);
    builder.add_dynamic(PySDM_dynamics.AmbientThermodynamics())
    builder.add_dynamic(PySDM_dynamics.Condensation())

    v_dry = formulae.trivia.volume(r_dry);
    r_wet = PySDM_initialisation.equilibrate_wet_radii(r_dry, environment, kappa * v_dry);
    attributes = py.dict(pyargs( ...
        'dry volume', formulae.trivia.volume(r_dry), ...
        'n', PySDM_initialisation.discretise_multiplicities(n_vol / rho0 * environment.mass_of_dry_air), ...
        'kappa times dry volume', kappa * v_dry, ...
        'volume', formulae.trivia.volume(r_wet) ...
    ));
    products = py.list({ ...
        PySDM_products.PeakSupersaturation(pyargs('name', 'S_max', 'unit', '%')), ...
        PySDM_products.AmbientDryAirDensity(), ...
        PySDM_products.ActivableFraction() ...
    });

    particulator = builder.build(pyargs( ...
        'attributes', attributes, ...
        'products', products ...
    ));
    P = particulator.products;
    get_value(P, 'S_max');
    S_max = -1;
    for i = 1:input.nt
        try
            particulator.run(int32(i));
        catch e
            e.message
            output = struct('N_act', nan, 'S_max', nan);
            return
        end
        S_max_curr = get_value(P, 'S_max');
        if S_max_curr < S_max
            break
        else
            S_max = S_max_curr;
        end
    end
    
    get = py.getattr(P{'activable fraction'}.get(pyargs('S_max', S_max)), '__getitem__');
    N1 = get(int32(0));
    
    output = struct(...
        'N_act', ...
        N1 * spectrum.norm_factor, ...
        'S_max', ...
        1 + S_max / 100 ...
    );
end

function [value] = get_value(products, name)
    get = py.getattr(products{name}.get(), '__getitem__');
    value = get(int32(0));
end
