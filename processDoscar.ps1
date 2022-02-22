$DistanceSource = @("2.12  3.12" -split " +")
$AtomType = 'VSe2'

enum MagTypeSource {
    # FM
    # AFM
    NM
}
foreach ($MagType in [MagTypeSource].GetEnumNames() ) {
    $FileName = $MagType.ToLower() + ".xml";
    foreach ($Distance in $DistanceSource) {
        # For VSe2_Primary
        # $scpXmlStr = "cluster:/home/yyzhao/TMDCs/clusterTransfer/primitiveCell/" + $AtomType + "/" + $Distance + "/" + $MagType.ToLower() + "/*.xml"
        # $scpDosStr = "cluster:/home/yyzhao/TMDCs/clusterTransfer/primitiveCell/" + $AtomType + "/" + $Distance + "/" + $MagType.ToLower() + "/dos/DOSCAR"

        # For VS2, VTe2_Primary
        # $scpXmlStr = "cluster:/home/yyzhao/TMDCs/primitiveCellCalc/mag/opt/" + $AtomType + "/change_D/" + $Distance + "/" + $MagType.ToLower() + "/*.xml"
        # $scpDosStr = "cluster:/home/yyzhao/TMDCs/primitiveCellCalc/mag/opt/" + $AtomType + "/change_D/" + $Distance + "/" + $MagType.ToLower() + "/dos/DOSCAR"

        # For VSe2 SuperCell FM Part1(2.42 - 3.12)
        # $FileName = $Distance.ToString() + ".xml";
        # $scpXmlStr = "cluster:/home/yyzhao/" + $AtomType + "/formal/tPhase/plane/paperFormal/1Tbilayer/" + $MagType.ToLower() + "/rotation/21.7867892982618/change_D_Ueff1/" + $Distance + "/*.xml";
        # $scpDosStr = "cluster:/home/yyzhao/" + $AtomType + "/formal/tPhase/plane/paperFormal/1Tbilayer/" + $MagType.ToLower() + "/rotation/21.7867892982618/change_D_Ueff1/" + $Distance + "/dos/DOSCAR";

        # For VSe2 SuperCell AFM Part1 (2.42-3.12)
        # $FileName = $Distance.ToString() + ".xml";
        # $scpXmlStr = "cluster:/home/yyzhao/VSe2/formal/tPhase/plane/paperFormal/1Tbilayer/afm/afm1/rotation/21.7867892982618/change_D_Ueff1/" + $Distance + "/*.xml"
        # $scpDosStr = "cluster:/home/yyzhao/VSe2/formal/tPhase/plane/paperFormal/1Tbilayer/afm/afm1/rotation/21.7867892982618/change_D_Ueff1/" + $Distance + "/dos/DOSCAR";

        # For VSe2 SuperCell FM/AFM part2 fm(1.62, 1.82-2.32) afm(2.12-2.32)
        # $FileName = $Distance.ToString() + ".xml";
        # $scpXmlStr = "cluster:/home/yyzhao/TMDCs/clusterTransfer/VSe2/" + $MagType.ToLower() + "/" + $Distance + "/*.xml"
        # $scpDosStr = "cluster:/home/yyzhao/TMDCs/clusterTransfer/VSe2/" + $MagType.ToLower() + "/" + $Distance + "/dos/DOSCAR";

        # For VS2 VTe2 Supercell
        # VTe2 "1.92  2.12  2.32  2.52  2.72  2.92  3.12  3.30  3.52"
        # $scpXmlStr = "cluster:/home/yyzhao/TMDCs/superCellCalc/" + $AtomType + "/" + $Distance + "/" + $MagType.ToLower() + "/*.xml"
        # $scpDosStr = "cluster:/home/yyzhao/TMDCs/superCellCalc/" + $AtomType + "/" + $Distance + "/" + $MagType.ToLower() + "/dos/DOSCAR";

        # # For VSe2 NM Primitivecell
        # # VSe2 "1.92  2.12  2.32  3.12"
        # $FileName = $Distance.ToString() + ".xml";
        # $scpXmlStr = "cluster:/home/yyzhao/TMDCs/nmCalc/VSe2/primitiveCell/" + $Distance + "/*.xml"
        # $scpDosStr = "cluster:/home/yyzhao/TMDCs/nmCalc/VSe2/primitiveCell/" + $Distance + "/dos/DOSCAR";

        # For VSe2 NM SuperCell
        # VSe2 "2.12  3.12"
        $FileName = $Distance.ToString() + ".xml";
        $scpXmlStr = "cluster:/home/yyzhao/TMDCs/nmCalc/VSe2/superCell/" + $Distance + "/*.xml"
        $scpDosStr = "cluster:/home/yyzhao/TMDCs/nmCalc/VSe2/superCell/" + $Distance + "/dos/DOSCAR";

        Write-Host $scpXmlStr;
        scp $scpXmlStr ./;
        [xml] $File = Get-Content $FileName;
        $DosPart = $File.modeling.summary | Where-Object {$_.name -match "dos"};
        # $EnergyPerAtom = $DosPart.totalEnergy/$DosPart.totAtomNum;
        # $TotalMagmom = $DosPart.totalmag;
        $TotalAtomNumber = $DosPart.totAtomNum;
        # $magmom = @($DosPart.magmont -split "\n");
        if ($TotalAtomNumber -eq 6) {
            $Struct = 'Primary'
        }
        else {
            $Struct = "Super"
        }

        # foreach ($cells in $magmom) {
        #     if ($cells -notmatch "tot") {
        #         @($AtomType, $Struct, $Distance, $MagType, $EnergyPerAtom, $TotalMagmom ,$cells) -join " " | Out-File -Append -FilePath .\databaseEng.dat -Encoding utf8
        #         # Write-Host $Struct, $Distance, $MagType, $EnergyPerAtom, $TotalMagmom ,$cells
        #     }    
        # }
        Write-Host $scpDosStr;
        scp $scpDosStr ./;
        .\processDoscar.exe $AtomType $Struct $Distance $MagType
    }
}

