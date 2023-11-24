{ lib, stdenv, fetchFromGitHub, cmake, gsl, gmp, makeWrapper, jemalloc, htslib, git, zlib, pkg-config }:

stdenv.mkDerivation rec {
  pname = "wfmash";
  version = "0.12.2";

  src = fetchFromGitHub {
    owner = "waveygang"; 
    repo = "wfmash";
    rev = "180455df2aaa298fbe2dafb1ddfc2796c78b3fcd";
    sha256 = "sha256-LxmVZcdOYWr1taMuzoW9/7bXJNo5NVgqOCz08gG7CuE=";
  };

  nativeBuildInputs = [ cmake makeWrapper ];

  buildInputs = [ 
    gsl 
    gmp
    jemalloc
    htslib
    git
    zlib
    pkg-config
  ];

  postPatch = ''
    mkdir -p include
    echo "#define WFMASH_GIT_VERSION \"${version}\"" > include/wfmash_git_version.hpp
  '';

  postInstall = ''
    wrapProgram $out/bin/wfmash --prefix PATH : ${lib.makeBinPath [ gsl gmp ]}
  '';

  meta = with lib; {
    description = "Base-accurate DNA sequence alignments using WFA and mashmap2";
    homepage = "https://github.com/ekg/wfmash";
    license = licenses.mit;
    platforms = platforms.linux;
    maintainers = [ maintainers.bzizou ];
  };
}
