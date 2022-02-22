$DistanceSource = @("1.92  2.12  2.32  2.52  2.72  2.92  3.12  3.30  3.52" -split " +")
$AtomType = 'VTe2'

enum MagTypeSource {
    FM
    AFM
}
foreach ($MagType in [MagTypeSource].GetEnumNames() ) {
    $FileName = $MagType.ToLower() + ".xml";
    foreach ($Distance in $DistanceSource) {
        # For VSe2_Primary
        # $scpStr = "cluster:/home/yyzhao/TMDCs/clusterTransfer/primitiveCell/" + $AtomType + "/" + $Distance + "/" + $MagType.ToLower() + "/*.xml"

        # For VS2, VTe2_Primary
        # $scpStr = "cluster:/home/yyzhao/TMDCs/primitiveCellCalc/mag/opt/" + $AtomType + "/change_D/" + $Distance + "/" + $MagType.ToLower() + "/*.xml"

        # For VSe2 SuperCell FM Part1
        # $FileName = $Distance.ToString() + ".xml";
        # $scpStr = "cluster:/home/yyzhao/" + $AtomType + "/formal/tPhase/plane/paperFormal/1Tbilayer/" + $MagType.ToLower() + "/rotation/21.7867892982618/change_D_Ueff1/" + $Distance + "/*.xml"

        # For VSe2 SuperCell AFM Part1
        # $FileName = $Distance.ToString() + ".xml";
        # $scpStr = "cluster:/home/yyzhao/VSe2/formal/tPhase/plane/paperFormal/1Tbilayer/afm/afm1/rotation/21.7867892982618/change_D_Ueff1/" + $Distance + "/*.xml"

        # For VSe2 SuperCell FM/AFM part2
        # $FileName = $Distance.ToString() + ".xml";
        # $scpStr = "cluster:/home/yyzhao/TMDCs/clusterTransfer/VSe2/" + $MagType.ToLower() + "/" + $Distance + "/*.xml"

        # For VS2 VTe2 Supercell
        $scpStr = "cluster:/home/yyzhao/TMDCs/superCellCalc/" + $AtomType + "/" + $Distance + "/" + $MagType.ToLower() + "/*.xml"
        Write-Host $scpStr;
        scp $scpStr ./;
        [xml] $File = Get-Content $FileName;
        $DosPart = $File.modeling.summary | Where-Object {$_.name -match "dos"};
        $EnergyPerAtom = $DosPart.totalEnergy/$DosPart.totAtomNum;
        $TotalMagmom = $DosPart.totalmag;
        $TotalAtomNumber = $DosPart.totAtomNum;
        $magmom = @($DosPart.magmont -split "\n");
        if ($TotalAtomNumber -eq 6) {
            $Struct = 'Primary'
        }
        else {
            $Struct = "Super"
        }

        foreach ($cells in $magmom) {
            if ($cells -notmatch "tot") {
                @($AtomType, $Struct, $Distance, $MagType, $EnergyPerAtom, $TotalMagmom ,$cells) -join " " | Out-File -Append -FilePath .\databaseEng.dat -Encoding utf8
                # Write-Host $Struct, $Distance, $MagType, $EnergyPerAtom, $TotalMagmom ,$cells
            }    
        }

    }
}

