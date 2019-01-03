sub generate_random_string{
    #http://stackoverflow.com/questions/13687643/generate-unique-random-strings
    my $length_of_randomstring = shift; # the length of 
                                        # the random string to generate

    my @chars=('a'..'z','A'..'Z','0'..'9','_');
    my $random_string;
    for(1..$length_of_randomstring){
        # rand @chars will generate a random 
        # number between 0 and scalar @chars
        $random_string.=$chars[rand @chars];
    }

    return $random_string;
}

