while (defined ($file = glob ("*.out"))) {
	$name = $file;
	$name =~ s/\..+//;
	$name = $name."_re.embl";
	open (In, "$file");
	open (Out, ">./$name");
		<In>;
		<In>;
		<In>;
		while (<In>) {
			if ($_ =~ /^\s/) {
				@temp = split (/\s+/, $_);
				$a = $temp[6]."..".$temp[7];
				$str = $temp[10].";".$temp[11];
				if ($temp[11] eq 'Low_complexity' || $temp[11] eq 'Simple_repeat') 	{next;}
				if ($temp[9] =~ /\+/) {
					$temp[14] =~ s/\D//g;
					$len = $temp[13] + $temp[14];
					$p = ($temp[13] - $temp[12] + 1) / $len;
					$note = "$temp[1];$temp[9];$temp[12];$temp[13];$len;$p";
						if ($temp[11]=~ /^LTR/) {
							print Out "FT   LTR             $a\n";
							print Out "FT                   /note=\"$str\"\n";
							print Out "FT                   /note=\"$note\"\n";
						}
						else {
							print Out "FT   repeat_region   $a\n";
							print Out "FT                   /note=\"$str\"\n";
							print Out "FT                   /note=\"$note\"\n";
						}
				}
				else {
					$temp[9] = "-";
					$temp[12] =~ s/\D//g;
					$len = $temp[13] + $temp[12];
					$p = ($temp[13] - $temp[14] + 1) / $len;
					$note = "$temp[1];$temp[9];$temp[14];$temp[13];$len;$p";
							if ($temp[11]=~ /^LTR/) {
							print Out "FT   LTR             complement($a)\n";
							print Out "FT                   /note=\"$str\"\n";
							print Out "FT                   /note=\"$note\"\n";
						}
							else {
							print Out "FT   repeat_region   complement($a)\n";
							print Out "FT                   /note=\"$str\"\n";
							print Out "FT                   /note=\"$note\"\n";
						}
				}

			}
			else {
				@temp = split (/\s+/, $_);
				$a = $temp[5]."..".$temp[6];
				$str = $temp[9].";".$temp[10];
				if ($temp[10] eq 'Low_complexity' || $temp[10] eq 'Simple_repeat') 	{next;}
				if ($temp[8] =~ /\+/) {
					$temp[13] =~ s/\D//g;
					$len = $temp[12] + $temp[13];
					$p = ($temp[12] - $temp[11] + 1) / $len;
					$note = "$temp[0];$temp[8];$temp[11];$temp[12];$len;$p";
							if ($temp[10]=~ /^LTR/) {
							print Out "FT   LTR             $a\n";
							print Out "FT                   /note=\"$str\"\n";
							print Out "FT                   /note=\"$note\"\n";
						}
							else {
							print Out "FT   repeat_region   $a\n";
							print Out "FT                   /note=\"$str\"\n";
							print Out "FT                   /note=\"$note\"\n";
						}
				}
				else {
					$temp[8] = "-";
					$temp[11] =~ s/\D//g;
					$len = $temp[12] + $temp[11];
					$p = ($temp[12] - $temp[13] + 1) / $len;
					$note = "$temp[0];$temp[8];$temp[13];$temp[12];$len;$p";
							if ($temp[10]=~ /^LTR/) {
							print Out "FT   LTR             complement($a)\n";
							print Out "FT                   /note=\"$str\"\n";
							print Out "FT                   /note=\"$note\"\n";
						}
							else {
							print Out "FT   repeat_region   complement($a)\n";
							print Out "FT                   /note=\"$str\"\n";
							print Out "FT                   /note=\"$note\"\n";
						}
				}
				
			}
		}
	close (Out);
	close (In);
}