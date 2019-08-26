n = 2; % number of images
res = 1024; % resolution
level = 7; % grid level
dir = 'classificationLetters/singleStepRP';

pix = 2^level+1;
for k = 0:n-2
    %% read in level set
    % fid = fopen( [dir, '/objectASCII', int2str( k ), '.pgm'] );
    % fseek( fid, 106, 'bof' );
    % levelset = fscanf( fid, '%g ' );
    % fclose( fid );
    levelset = double( imread( [dir, '/object0_', int2str( k ), '.pgm'], 'PGM' ) );
    % [dir, '/object', int2str( k ), '.pgm'],
    % levelset = textread( [dir, '/object0_', int2str( k ), '.pgm'], '', 'headerlines', 4, 'whitespace', ' \b\t\n' );
    % levelset = levelset( 1: size( levelset, 2 ) - 1 );
    % levelset = reshape( levelset, 129, 129 )';
    display('read in object');

    %% read in energy distribution
    % fid = fopen( [dir, '/elasticEnergyASCII', int2str( k ), '.pgm'] );
    % fseek( fid, 111, 'bof' );
    % energy = fscanf( fid, '%g ' );
    % fclose( fid );
    energy = textread( [dir, '/elastEnergy', int2str( k ), '.pgm'], '', 'headerlines', 4, 'whitespace', ' \b\t\n' );
    energy = energy( 1: size( energy, 2 ) - 1 );
    energy = reshape( energy, pix, pix )';
    display('read in energy');

    %% resample files
    x = 1:pix;
    y = x;
    [x,y] = meshgrid( x, y );
    X = linspace( 1, pix, res );
    Y = X;
    [X,Y] = meshgrid( X, Y );
    Levelset = interp2( x, y, levelset, X, Y );
    Energy = interp2( x, y, energy, X, Y );
    % surf( X, Y, Energy, min( Energy, .3 ) );
    % view(2);
    % shading interp;
    % axis equal;
    % axis( [ 1, pix, 1, pix ] );
    display('interpolated object and energy');

    %% save as png
    Levelset = Levelset - min( min( Levelset ) );
    Levelset = Levelset / max( max( Levelset ) );
    Levelset = ( .5 < Levelset ) * 1.;
    jet = colormap( jet );
    %energyVals = linspace( 0, 3., size( jet, 1 ) );
    %energyVals = linspace( 0, .03, size( jet, 1 ) );
    energyVals = logspace( .5, 2, size( jet, 1 ) )-10^.5;
    energyIm = zeros( size( Energy, 1 ), size( Energy, 2 ), 3 );
    for i = 1:3
        energyIm(:,:,i) = interp1( energyVals, jet(:,i), Energy, 'linear', jet( size( jet, 1 ), i ) );
    end
    imwrite( energyIm, [dir, '/object', int2str( k ), '.png'], 'png', 'Alpha', Levelset );

end