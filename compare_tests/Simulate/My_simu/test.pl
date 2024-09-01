#! /usr/bin/perl -w
my $repeat_var = "a";
my $transition = "b";
my $repeat_var = "c";
my $var = $repeat_var.$transition.$repeat_var;
ins_seq("ASDAD");
sub insert_seq {
    my $var = shift;
    say $var
}