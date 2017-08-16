
function [restrict, indeks] = restrykcje(n, t, vectors, vec_obs, b_C, war, dim)

switch war.case
    case 1
        restrict = zeros(size(vectors,2),2,dim(2)); %bo  po dwie restrykcje
        %w zmiennej indeks kolumnami poustawiane sa wszystkie zmienne
        %po 2 restrykcje na zmienna
        %dim(2) czyli tyle restrykcji ile zmiennych
        indeks = zeros(2,dim(2));
        for a = 1:dim(2)
            %disp(size(vectors));
            elast = vectors(:,:,a)*b_C;
            [value, index] = min(elast);
            indeks(1,a) = index;
            fprintf('Dla %d zmiennej znaleziono 1 minimum: %d \n',a,value);
            restrict(:,1,a) = vectors(index,:,a)';
            [value, index2] = max(elast);
            fprintf('Dla %d zmiennej znaleziono max: %d na podmiane\n',a,value);
            %po co miala byc ta podmiana?
            %odp: zeby po raz drugi nie wykryl mi tego samego pkt
            elast(index) = elast(index2);
            [value, index] = min(elast);
            indeks(2,a) = index;
            fprintf('Dla %d zmiennej znaleziono 2 minimum: %d \n',a,value);
            restrict(:,2,a) = vectors(index,:,a)';
        end
        restrict = restrict(:,:)';
        indeks = indeks(:);
    case 2 %po dwie restrykcje na zmienna dla 'm' pierwszych objektow
        zmien = size(vectors);
        m = war.m;
        diff = n-m;
        vectors_new = zeros(m*t,zmien(2),zmien(3));
        for b = 1:zmien(3)
            a1 = 1;
            a11 = 1;
            for a2 = 1:t
                for a3 = 1:m
                    %disp(size(restrict(a1,:)));
                    %disp(size(vectors(a11,:,war.zm)'));
                    vectors_new(a1,:,b) = vectors(a11,:,b);
                    a1 = a1+1;
                    a11 = a11+1;
                end
                a11 = a11+diff;
            end
        end
        
        restrict = zeros(size(vectors_new,2),2,dim(2)); 
        indeks = zeros(2,dim(2));
        for a = 1:dim(2)
            %disp(size(vectors));
            elast = vectors_new(:,:,a)*b_C;
            [value, index] = min(elast);
            indeks(1,a) = index;
            fprintf('Dla %d zmiennej znaleziono 1 minimum: %d \n',a,value);
            restrict(:,1,a) = vectors_new(index,:,a)';
            [value, index2] = max(elast);
            fprintf('Dla %d zmiennej znaleziono max: %d na podmiane\n',a,value);
            %po co miala byc ta podmiana?
            %odp: zeby po raz drugi nie wykryl mi tego samego pkt
            elast(index) = elast(index2);
            [value, index] = min(elast);
            indeks(2,a) = index;
            fprintf('Dla %d zmiennej znaleziono 2 minimum: %d \n',a,value);
            restrict(:,2,a) = vectors_new(index,:,a)';
        end
        restrict = restrict(:,:)';
        indeks = indeks(:);
        
        
    case 3 %restrykcja na jedna ze zmiennych (w calosci)
        %wymagany war.zm do zdefiniowania
        restrict = vectors(:,:,war.zm);
        indeks = 0;
    case 4 %restrykcje na jedna zmienna dla m pierwszych objektow
        m = war.m;
        restrict = zeros(m*t,size(vectors(1,:,1),2));
        diff = n-m;
        a1 = 1;
        a11 = 1;
        for a2 = 1:t
            for a3 = 1:m
                %disp(size(restrict(a1,:)));
                %disp(size(vectors(a11,:,war.zm)'));
                restrict(a1,:) = vectors(a11,:,war.zm)';
                a1 = a1+1;
                a11 = a11+1;
            end
            a11 = a11+diff;
        end
        indeks =0;
        
    case 5 %restrykcje dla srednich po obserwacjach, po jednej restrykcji
        restrict = zeros(size(vec_obs,2),1,dim(2)); %bo  po jednej restrykcji
        %w zmiennej indeks kolumnami poustawiane sa wszystkie zmienne
        %po 2 restrykcje na zmienna
        %dim(2) czyli tyle restrykcji ile zmiennych
        indeks = zeros(1,dim(2));
        for a = 1:dim(2)
            elast = vec_obs(:,:,a)*b_C;
            [value, index] = min(elast);
            indeks(1,a) = index;
            fprintf('Dla %d zmiennej znaleziono 1 minimum: %d \n',a,value);
            restrict(:,1,a) = vec_obs(index,:,a)';
        end
        restrict = restrict(:,:)';
        indeks = indeks(:);
    case 6 %restrykcje dla srednich po obserwacjach, po jednej restrykcji
        %dla pierwszych m objektow
        vec_obs = vec_obs(1:war.m,:,:);
        restrict = zeros(size(vec_obs,2),1,dim(2)); %bo  po jednej restrykcji
        %w zmiennej indeks kolumnami poustawiane sa wszystkie zmienne
        %po 2 restrykcje na zmienna
        %dim(2) czyli tyle restrykcji ile zmiennych
        indeks = zeros(1,dim(2));
        for a = 1:dim(2)
            elast = vec_obs(:,:,a)*b_C;
            [value, index] = min(elast);
            indeks(1,a) = index;
            fprintf('Dla %d zmiennej znaleziono 1 minimum: %d \n',a,value);
            restrict(:,1,a) = vec_obs(index,:,a)';
        end
        restrict = restrict(:,:)';
        indeks = indeks(:);
    case 7 %restrykcje dla srednich po obserwacjach, wszystkich
        restrict = zeros(size(vec_obs,2),n,dim(2)); %bo  po jednej restrykcji
        %dim(2) czyli tyle restrykcji ile zmiennych
        for a = 1:dim(2)
            %disp(size(restrict(:,:,a)));
            %disp(size(vec_obs(:,:,a)'));
            restrict(:,:,a) = vec_obs(:,:,a)';
        end
        restrict = restrict(:,:)';
        indeks = 0;
    case 8 %restrykcje dla srednich po obserwacjach, wszystkich
        %dla pierwszych m objektow
        vec_obs = vec_obs(1:war.m,:,:);
        restrict = zeros(size(vec_obs,2),war.m,dim(2)); %bo  po jednej restrykcji
        %dim(2) czyli tyle restrykcji ile zmiennych
        for a = 1:dim(2)
            %disp(size(restrict(:,:,a)));
            %disp(size(vec_obs(:,:,a)'));
            restrict(:,:,a) = vec_obs(:,:,a)';
        end
        restrict = restrict(:,:)';
        indeks = 0;
        
    otherwise
        disp('nie ma takiego warunku szwagier');
end