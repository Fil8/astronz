# AstronZ

**Insallation instructions**

In a new virtual environment, clone this repository: `git clone git clone https://github.com/Fil8/astronz.
From the directory of your repository run:

```
pip install -e .
```

***

**Requisites**

AstronZ provides a set of tools for the analysis of astronomical observations. Written in python, it makes use of the most common libraries (e.g. `numpy`, `scipy`, `astropy`). 


***
**Description**

This is a set of tools that can be useful for quick calculations of astronomical parameters. Three main classes of functions are available:


- **Cosmo** : set of cosmological tools
- **radioHI** : set of tools useful for neutral hydrogen (_HI_) science
- **AGN** : set of tools useful for general AGN science

check out the [WiKi](https://github.com/Fil8/astronz/wiki) for a complete description of all functions available in the `astronz` classes.

***
**Usage**

AstronZ takes its variables from a default parameter file and from terminal, if any are given. 

From your current working directory typying `astronz` this message will be shown: 

```
	************* --- AstronZ --- **************

	   ... list of the avaliable classes: ...

	 - c (cosmological tools)
	 - hi (neutral hydrogen tools)
	 - a (AGN science tools)


```

Typing one of the available commands (`c`, `hi` or `a`) you will be redirected to the chosen class where you can select the function of interest. Ingesting from terminal the required parameters, `astronz` will output the result on screen. 

This is an example, where using the class `hi` we determine the frequency of the neutral hydrogen line at redshift, $z=1$:


```
(astronz) maccagni@monaco:~/programs/astronz$ astronz

	************* --- AstronZ --- **************

	   ... list of the avaliable classes: ...

	 - c (cosmological tools)
	 - hi (neutral hydrogen tools)
	 - a (AGN science tools)

>>> hi

	************* --- AstronZ : HI --- **************

		... Neutral Hydrogen Tools ... 

	... Here's what I can do: ...

	- HI at z	->	HI
	- Tully-Fisher	->	TF
	- TF-apparent	->	TFM
	- Absolute Magnitude	->	M
	- Optical depth	->	TAU
	- Abs. Column Density	->	NHI
	- fluxColumn Density	->	FHI
	- EMColumn Density	->	EHI
	- Mass HI	->	MHI
	- Mass HI	->	CHI
	- Mimimum_mass	->	MM
	- Radio Power	->	RP
	- Vel. at z(HI)	->	VEL
	- Incl. deVAB	->deVAB
	- Vel. res.	->	VRES
        
>>> hi

z= 1
HI = 710.203 MHz

	************* --- HI : DONE --- **************

	************ --- AstronZ : DONE --- *************
```

**Other usages**

You can directly access `astronz` classes by typing directly one of the following commands: 

- `astronz -c` : for cosmological tools
- `astronz -hi`: for neutral hydrogen tools
- `astronz -a` : for AGN science tools

**Help**

`astronz -h` will show you a (minimal) help:

```
	************* --- AstronZ : Help --- **************

		  ... called for help ...

usage: astronz [-h] [-v] [-c] [-a] [-hi]

AstronZ: tools to analyse astronomical data

version 1.0.0

install path /home/maccagni/programs/astronz/astronz

Filippo Maccagni <filippo.maccagni@gmial.com>

optional arguments:
  -h, --help      Print help message and exit
  -v, --version   show program's version number and exit
  -c, --cosmo     tools for cosmological calculations
  -a, --agn       tools for AGN science
  -hi, --radioHI  tools for neutral hydrogen science

Run a command. This can be:

astronz		(all tools)
astronz -c	(cosmological tools)
astronz -hi	(neutral hydroge tools)
astronz -agn 	(AGN science tools)
            

	************* --- AstronZ : DONE --- **************

```

***

**License**

This project is licensed under the GNU General Public License v3.0 - see [license](https://github.com/Fil8/astronz/blob/master/LICENSE.md) for details.


 ***
 <p>&copy <sub> Filippo M. Maccagni 2018 </sub></p>
