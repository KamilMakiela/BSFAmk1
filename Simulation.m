function [wyniki, pliki] = Simulation(y, x, n, t, u, beta, iteracji, model, H, objetosc, spalone, restrykcje, rozklad, mediana, priors)

switch rozklad
    case 1
    switch model
        case 1
            [wyniki, pliki] = SimulationER1Simple(y, x, n, t, u, beta, iteracji, H, objetosc, spalone, restrykcje, mediana, priors);
        case 2
            [wyniki, pliki] = SimulationER1Simple(y, x, n, t, u, beta, iteracji, H, objetosc, spalone, restrykcje, mediana, priors);
        case 4
            [wyniki, pliki] = SimulationER1Simple(y, x, n, t, u, beta, iteracji, H, objetosc, spalone, restrykcje, mediana, priors);
        case {6, 8, 10, 12}
            restrykcje.restrict = NaN;
            [wyniki, pliki] = SimulationER1Simple(y, x, n, t, u, beta, iteracji, H, objetosc, spalone, restrykcje, mediana, priors);
        case {5,7, 9, 11}
            [wyniki, pliki] = SimulationER1Alt(y, x, n, t, u, beta, iteracji, H, objetosc, spalone, restrykcje, mediana, priors);
        otherwise
            msgbox('Unknown model','Well this is embarrassing...','error');
            wyniki = 'Ups';
    end
    case 2
    switch model
        case 1
            [wyniki, pliki] = SimulationHNSimple(y, x, n, t, u, beta, iteracji, H, objetosc, spalone, restrykcje, mediana, priors);
        case 2
            [wyniki, pliki] = SimulationHNSimple(y, x, n, t, u, beta, iteracji, H, objetosc, spalone, restrykcje, mediana, priors);
        case 4
            [wyniki, pliki] = SimulationHNSimple(y, x, n, t, u, beta, iteracji, H, objetosc, spalone, restrykcje, mediana, priors);
        case {6, 8, 10, 12}
            restrykcje.restrict = NaN;
            [wyniki, pliki] = SimulationHNSimple(y, x, n, t, u, beta, iteracji, H, objetosc, spalone, restrykcje, mediana, priors);
        case {5,7, 9, 11}
            [wyniki, pliki] = SimulationHNAlt(y, x, n, t, u, beta, iteracji, H, objetosc, spalone, restrykcje, mediana, priors);
        otherwise
            msgbox('Unknown model','Well this is embarrassing...','error');
            wyniki = 'Ups';
    end
    case {3,4}
    switch model
        case 1
            [wyniki, pliki] = SimulationVED(y, x, n, t, u, beta, iteracji, H, objetosc, spalone, restrykcje, mediana, priors);
        case 2
            [wyniki, pliki] = SimulationVED(y, x, n, t, u, beta, iteracji, H, objetosc, spalone, restrykcje, mediana, priors);
        case 4
            [wyniki, pliki] = SimulationVED(y, x, n, t, u, beta, iteracji, H, objetosc, spalone, restrykcje, mediana, priors);
        case {6, 8, 10, 12}
            restrykcje.restrict = NaN;
            [wyniki, pliki] = SimulationVED(y, x, n, t, u, beta, iteracji, H, objetosc, spalone, restrykcje, mediana, priors);
        case {5,7, 9, 11}
            [wyniki, pliki] = SimulationVED(y, x, n, t, u, beta, iteracji, H, objetosc, spalone, restrykcje, mediana, priors);
        otherwise
            msgbox('Unknown model','Well this is embarrassing...','error');
            wyniki = 'Ups';
    end
    case {5}
        [wyniki, pliki] = SimulationProsty(y, x, n, t, u, beta, iteracji, H, objetosc, spalone, restrykcje, priors);
    
    otherwise 
        msgbox('Unknown model.');
        wyniki = 'Ups';
end
